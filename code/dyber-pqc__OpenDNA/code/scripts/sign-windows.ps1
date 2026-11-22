<#
.SYNOPSIS
    Sign one Windows binary during a Tauri build.

.DESCRIPTION
    Tauri calls this once per executable and installer it produces, passing the
    file path. It is wired in by app/desktop/tauri.signing.conf.json, which is
    merged into the build only when a signing credential is present, so an
    ordinary unsigned build never reaches this script.

    There are two ways to tell it what to do, and which you need depends on
    where the private key lives.

    WINDOWS_CERT_THUMBPRINT
        The SHA-1 thumbprint of a certificate already in the machine's
        certificate store, signed with signtool.exe from the Windows SDK. This
        is the path for a certificate on a hardware token in front of a person,
        or an old file-based certificate. It does not suit CI, because the key
        has to be in the store on the runner.

    WINDOWS_SIGN_COMMAND
        A command template with {} standing for the file to sign. This is the
        path for a cloud signing service, which is what a new publicly-trusted
        certificate now requires: since June 2023 the private key has to live
        on certified hardware, so you rent signing from a service rather than
        holding a .pfx. A template rather than built-in support for one vendor,
        because the vendors' tools differ and this script should not pretend to
        know flags it has never run. docs/SIGNING.md has worked examples.

    With neither set the script fails loudly. A build that silently produced
    unsigned binaries while claiming to sign them would be worse than one that
    stops.

.PARAMETER Path
    The file to sign, passed by Tauri.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Path
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $Path)) {
    throw "sign-windows: no such file: $Path"
}

# RFC 3161 timestamping, so signatures stay valid after the certificate
# expires. Overridable because a vendor may require its own.
$timestampUrl = if ($env:WINDOWS_TIMESTAMP_URL) { $env:WINDOWS_TIMESTAMP_URL } else { 'http://timestamp.digicert.com' }

if ($env:WINDOWS_SIGN_COMMAND) {
    # The template owns its own quoting; {} is replaced with the path.
    $command = $env:WINDOWS_SIGN_COMMAND.Replace('{}', $Path)
    Write-Host "sign-windows: running the configured signing command on $Path"
    cmd.exe /c $command
    if ($LASTEXITCODE -ne 0) {
        throw "sign-windows: the signing command exited $LASTEXITCODE"
    }
}
elseif ($env:WINDOWS_CERT_THUMBPRINT) {
    $signtool = (Get-Command signtool.exe -ErrorAction SilentlyContinue).Source
    if (-not $signtool) {
        # The SDK installs it outside PATH; look where it actually lands.
        $candidates = Get-ChildItem -Path 'C:\Program Files (x86)\Windows Kits\10\bin' -Filter 'signtool.exe' -Recurse -ErrorAction SilentlyContinue |
            Where-Object { $_.FullName -like '*\x64\*' } |
            Sort-Object FullName -Descending
        if ($candidates) { $signtool = $candidates[0].FullName }
    }
    if (-not $signtool) {
        throw 'sign-windows: signtool.exe not found. Install the Windows SDK or use WINDOWS_SIGN_COMMAND.'
    }
    Write-Host "sign-windows: signing $Path with thumbprint $($env:WINDOWS_CERT_THUMBPRINT)"
    & $signtool sign /fd sha256 /td sha256 /tr $timestampUrl /sha1 $env:WINDOWS_CERT_THUMBPRINT $Path
    if ($LASTEXITCODE -ne 0) {
        throw "sign-windows: signtool exited $LASTEXITCODE"
    }
}
else {
    throw @'
sign-windows: no signing credential configured.

Set one of:
  WINDOWS_SIGN_COMMAND   a command template, {} standing for the file to sign
  WINDOWS_CERT_THUMBPRINT  a certificate in the machine store, signed with signtool

This script is only reached when tauri.signing.conf.json has been merged into
the build, which the release workflow does only when a credential is present.
Reaching it without one means the workflow and the secrets disagree.
See docs/SIGNING.md.
'@
}

Write-Host "sign-windows: signed $Path"
