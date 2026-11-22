use serde::{Deserialize, Serialize};
use sysinfo::System;

/// Detected hardware capabilities.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HardwareInfo {
    /// CPU model name.
    pub cpu_name: String,
    /// Logical CPU count.
    pub cpu_cores: usize,
    /// Installed RAM in GiB.
    pub total_ram_gb: f64,
    /// Currently available RAM in GiB.
    pub available_ram_gb: f64,
    /// Primary GPU, if one was detected.
    pub gpu: Option<GpuInfo>,
    /// Operating system (`std::env::consts::OS`).
    pub os: String,
    /// Model tier this machine can run comfortably.
    pub recommended_tier: ModelTier,
    /// Best available compute backend.
    pub recommended_backend: ComputeBackend,
}

/// A detected GPU.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GpuInfo {
    /// GPU model name.
    pub name: String,
    /// Dedicated memory in GiB.
    pub vram_gb: f64,
    /// Backend used to drive this GPU.
    pub backend: ComputeBackend,
}

/// Compute backend for ML inference.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum ComputeBackend {
    /// NVIDIA CUDA.
    Cuda,
    /// Apple Metal.
    Metal,
    /// AMD `ROCm`.
    Rocm,
    /// Windows `DirectML`.
    DirectML,
    /// CPU only.
    Cpu,
}

impl std::fmt::Display for ComputeBackend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ComputeBackend::Cuda => write!(f, "CUDA (NVIDIA)"),
            ComputeBackend::Metal => write!(f, "Metal (Apple)"),
            ComputeBackend::Rocm => write!(f, "ROCm (AMD)"),
            ComputeBackend::DirectML => write!(f, "DirectML (Windows)"),
            ComputeBackend::Cpu => write!(f, "CPU"),
        }
    }
}

/// Model complexity tier based on hardware.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, PartialOrd)]
pub enum ModelTier {
    /// 4GB RAM, CPU only - smallest models
    Minimal = 1,
    /// 8GB VRAM, RTX 3060 - standard models
    Standard = 2,
    /// 16GB+ VRAM, RTX 4080/4090 - large models
    Large = 3,
    /// 24GB+ VRAM, A6000/H100 - full suite
    Full = 4,
}

impl std::fmt::Display for ModelTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ModelTier::Minimal => write!(f, "Minimal (CPU only, smallest models)"),
            ModelTier::Standard => write!(f, "Standard (GPU, standard models)"),
            ModelTier::Large => write!(f, "Large (GPU, large models)"),
            ModelTier::Full => write!(f, "Full (GPU, full model suite)"),
        }
    }
}

/// Inference limits derived from the hardware tier.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DeviceCapability {
    /// Longest sequence to fold in one pass.
    pub max_sequence_length: usize,
    /// Numeric precision for inference.
    pub recommended_precision: Precision,
    /// Sequences per batch.
    pub max_batch_size: usize,
    /// Whether fused attention kernels are available.
    pub supports_flash_attention: bool,
}

/// Numeric precision for model weights and activations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Precision {
    /// 32-bit float.
    Float32,
    /// 16-bit float.
    Float16,
    /// bfloat16.
    BFloat16,
    /// 8-bit integer quantisation.
    Int8,
}

impl std::fmt::Display for Precision {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Precision::Float32 => write!(f, "FP32"),
            Precision::Float16 => write!(f, "FP16"),
            Precision::BFloat16 => write!(f, "BF16"),
            Precision::Int8 => write!(f, "INT8"),
        }
    }
}

/// Detect hardware and recommend configuration.
///
/// # Errors
/// Returns the `sysinfo` error if system information cannot be read.
pub fn detect_hardware() -> Result<HardwareInfo, sysinfo::Error> {
    let sys = System::new_all()?;

    let cpu_name = sys.cpus().first().map_or_else(|| "Unknown CPU".to_string(), |c| c.brand().to_string());

    let cpu_cores = sys.cpus().len();
    let total_ram_gb = sys.total_memory() as f64 / 1_073_741_824.0;
    let available_ram_gb = sys.available_memory() as f64 / 1_073_741_824.0;

    let os = std::env::consts::OS.to_string();

    // GPU detection: nvidia-smi for CUDA, otherwise assume Apple Silicon on macOS.
    // ROCm and DirectML probing are not implemented yet.
    let gpu = detect_gpu(&os);

    let recommended_tier = determine_tier(gpu.as_ref(), total_ram_gb);
    let recommended_backend = gpu.as_ref().map_or(ComputeBackend::Cpu, |g| g.backend);

    Ok(HardwareInfo {
        cpu_name,
        cpu_cores,
        total_ram_gb,
        available_ram_gb,
        gpu,
        os,
        recommended_tier,
        recommended_backend,
    })
}

fn detect_gpu(os: &str) -> Option<GpuInfo> {
    // Try NVIDIA first (most common for ML)
    if let Some(gpu) = detect_nvidia() {
        return Some(gpu);
    }

    // On macOS, assume Apple Silicon with unified memory
    if os == "macos" {
        return Some(GpuInfo {
            name: "Apple Silicon (unified memory)".to_string(),
            vram_gb: 0.0, // Will use system RAM
            backend: ComputeBackend::Metal,
        });
    }

    None
}

fn detect_nvidia() -> Option<GpuInfo> {
    // Try nvidia-smi to detect GPU
    let output = std::process::Command::new("nvidia-smi")
        .args(["--query-gpu=name,memory.total", "--format=csv,noheader,nounits"])
        .output()
        .ok()?;

    if !output.status.success() {
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let line = stdout.lines().next()?;
    let parts: Vec<&str> = line.split(',').map(str::trim).collect();

    if parts.len() >= 2 {
        let name = parts[0].to_string();
        let vram_mb: f64 = parts[1].parse().unwrap_or(0.0);
        Some(GpuInfo { name, vram_gb: vram_mb / 1024.0, backend: ComputeBackend::Cuda })
    } else {
        None
    }
}

fn determine_tier(gpu: Option<&GpuInfo>, total_ram_gb: f64) -> ModelTier {
    match gpu {
        Some(g) if g.vram_gb >= 24.0 => ModelTier::Full,
        Some(g) if g.vram_gb >= 16.0 => ModelTier::Large,
        Some(g) if g.vram_gb >= 6.0 => ModelTier::Standard,
        Some(g) if g.backend == ComputeBackend::Metal && total_ram_gb >= 16.0 => ModelTier::Standard,
        Some(g) if g.backend == ComputeBackend::Metal => ModelTier::Minimal,
        _ => ModelTier::Minimal,
    }
}

/// Inference limits for the detected hardware.
pub fn get_device_capability(hw: &HardwareInfo) -> DeviceCapability {
    match hw.recommended_tier {
        ModelTier::Minimal => DeviceCapability {
            max_sequence_length: 256,
            recommended_precision: Precision::Float32,
            max_batch_size: 1,
            supports_flash_attention: false,
        },
        ModelTier::Standard => DeviceCapability {
            max_sequence_length: 512,
            recommended_precision: Precision::Float16,
            max_batch_size: 4,
            supports_flash_attention: true,
        },
        ModelTier::Large => DeviceCapability {
            max_sequence_length: 1024,
            recommended_precision: Precision::Float16,
            max_batch_size: 8,
            supports_flash_attention: true,
        },
        ModelTier::Full => DeviceCapability {
            max_sequence_length: 2048,
            recommended_precision: Precision::Float16,
            max_batch_size: 16,
            supports_flash_attention: true,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_hardware() {
        let hw = detect_hardware().unwrap();
        assert!(hw.cpu_cores > 0);
        assert!(hw.total_ram_gb > 0.0);
    }

    #[test]
    fn test_tier_determination() {
        let gpu_24gb =
            Some(GpuInfo { name: "RTX 4090".to_string(), vram_gb: 24.0, backend: ComputeBackend::Cuda });
        assert_eq!(determine_tier(gpu_24gb.as_ref(), 32.0), ModelTier::Full);

        assert_eq!(determine_tier(None, 8.0), ModelTier::Minimal);
    }
}
