<#
.SYNOPSIS
    Run cargo with the MSVC and CUDA toolchains on PATH, for GPU builds.

.DESCRIPTION
    Building candle's CUDA backend on Windows needs two things that are not
    on PATH by default: nvcc from the CUDA toolkit, and cl.exe from Visual
    Studio's C++ tools (nvcc shells out to it). This script loads both, then
    forwards its arguments to cargo.

    candle 0.11 supports CUDA up to 13.4; newer toolkits fail the build with
    "Unsupported cuda toolkit version".

    Do not redirect this script's output (`2>$null`, `*>file`, `| tee`).
    Cargo writes its progress to stderr, and PowerShell turns a native
    program's stderr into NativeCommandError records, so the run reports
    failure even when cargo exited 0. Let it write to the console, or capture
    it from outside PowerShell.

    Prefer long flags for the program being run (`--input`, not `-i`). This
    script has CmdletBinding, so PowerShell tries to bind a dash argument to
    its own common parameters before forwarding it, and anything that is a
    unique prefix of one is swallowed: `-i` fails with "the parameter name
    'i' is ambiguous" against -InformationAction and -InformationVariable.
    Letters matching no common parameter, such as cargo's `-j`, pass through
    unharmed, which is why the first example below works.

    Quote the `--` separator as '--'. Bare, PowerShell consumes it as its own
    end-of-parameters token and cargo never sees it, so the program's first
    flag is read as one of cargo's and rejected.

    Compiling the CUDA kernels starts one nvcc per core and each is a heavy
    C++ compile; on a machine with many cores and limited RAM that can
    exhaust memory and lock up the desktop. Pass `-j 4` (or lower) for the
    first build, after which the kernels are cached.

.EXAMPLE
    ./scripts/cargo-cuda.ps1 build -j 4 -p opendna-ml --features cuda --release

.EXAMPLE
    ./scripts/cargo-cuda.ps1 run -p opendna-train --features cuda '--' embed --task disorder
#>
[CmdletBinding()]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]] $CargoArgs
)

$ErrorActionPreference = 'Stop'

# The toolkit must satisfy two upper bounds:
#   * the driver's, because nvcc emits PTX that the driver has to load, and a
#     newer toolkit fails at runtime with CUDA_ERROR_UNSUPPORTED_PTX_VERSION;
#   * candle 0.11's, via cudarc 0.19, which rejects anything above 13.4 at
#     build time.
# Set CUDA_PATH yourself to override this choice.
$candleMax = [version]'13.4'

function Get-DriverCudaVersion {
    $smi = Get-Command nvidia-smi -ErrorAction SilentlyContinue
    if (-not $smi) { return $null }
    $line = (& $smi.Source 2>$null | Select-String -Pattern 'CUDA Version:\s*([0-9]+\.[0-9]+)')
    if ($line) { return [version]$line.Matches[0].Groups[1].Value }
    return $null
}

$driverMax = Get-DriverCudaVersion

# An inherited CUDA_PATH is often whichever toolkit was installed last, which
# is not necessarily one this driver can run. Check it before trusting it.
if ($env:CUDA_PATH -and $env:CUDA_PATH -match 'v([0-9]+\.[0-9]+)$') {
    $current = [version]$Matches[1]
    if ($current -gt $candleMax -or ($driverMax -and $current -gt $driverMax)) {
        Write-Warning "CUDA_PATH points at $current, which this setup cannot use; looking for another toolkit."
        $env:CUDA_PATH = $null
    }
}

if (-not $env:CUDA_PATH) {
    $root = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA'
    $installed = @()
    if (Test-Path $root) {
        $installed = Get-ChildItem $root -Directory |
            Where-Object { $_.Name -match '^v([0-9]+\.[0-9]+)$' } |
            ForEach-Object { [pscustomobject]@{ Version = [version]$Matches[1]; Path = $_.FullName } } |
            Sort-Object Version -Descending
    }
    $usable = $installed | Where-Object {
        $_.Version -le $candleMax -and (-not $driverMax -or $_.Version -le $driverMax)
    }
    if (-not $usable) {
        $have = if ($installed) { ($installed.Version -join ', ') } else { 'none' }
        $limit = if ($driverMax) { "driver supports up to $driverMax" } else { 'driver version unknown' }
        throw "No usable CUDA toolkit (installed: $have; $limit; candle supports up to $candleMax). " +
              "Install a matching one, e.g. winget install --id Nvidia.CUDA --version 13.2"
    }
    $env:CUDA_PATH = $usable[0].Path
}
if (-not (Test-Path $env:CUDA_PATH)) {
    throw "CUDA_PATH points at $env:CUDA_PATH, which does not exist."
}
# CUDA 13 moved the runtime DLLs (cudart64_13.dll, cublas64_13.dll, ...) into
# an architecture subdirectory; 12.x kept them directly in bin. Add both, or
# the compiled binary starts and immediately fails with STATUS_DLL_NOT_FOUND.
$env:PATH = "$env:CUDA_PATH\bin;$env:PATH"
$cudaBinX64 = Join-Path $env:CUDA_PATH 'bin\x64'
if (Test-Path $cudaBinX64) {
    $env:PATH = "$cudaBinX64;$env:PATH"
}

# Locate the C++ compiler through vswhere rather than hard-coding a version.
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) {
    throw "vswhere.exe not found; install Visual Studio with the C++ build tools."
}
$vsPath = & $vswhere -latest -products * `
    -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsPath) {
    throw "No Visual Studio installation with the C++ build tools was found."
}
$vcvars = Join-Path $vsPath 'VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) {
    throw "vcvars64.bat not found under $vsPath."
}

# vcvars64.bat only sets variables in its own shell, so run it and import
# the resulting environment.
& cmd.exe /c "`"$vcvars`" >nul 2>&1 && set" | ForEach-Object {
    if ($_ -match '^([^=]+)=(.*)$') {
        Set-Item -Path "env:$($Matches[1])" -Value $Matches[2] -ErrorAction SilentlyContinue
    }
}

if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
    throw "cl.exe is still not on PATH after loading $vcvars."
}

# CUDA 13's bundled CCCL headers refuse to compile under MSVC's traditional
# preprocessor. nvcc reads NVCC_APPEND_FLAGS itself, which is the only way to
# reach cl.exe here: candle-kernels builds its nvcc command line internally
# and exposes no hook for extra flags.
$zc = '-Xcompiler=/Zc:preprocessor'
if ($env:NVCC_APPEND_FLAGS -and $env:NVCC_APPEND_FLAGS -notlike "*$zc*") {
    $env:NVCC_APPEND_FLAGS = "$env:NVCC_APPEND_FLAGS $zc"
} elseif (-not $env:NVCC_APPEND_FLAGS) {
    $env:NVCC_APPEND_FLAGS = $zc
}

Write-Host "CUDA $(Split-Path $env:CUDA_PATH -Leaf), MSVC $((Get-Command cl.exe).Source)" -ForegroundColor DarkGray
& cargo @CargoArgs
exit $LASTEXITCODE
