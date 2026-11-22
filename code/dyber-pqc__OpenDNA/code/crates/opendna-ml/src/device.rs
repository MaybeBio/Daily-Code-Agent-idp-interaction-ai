use candle_core::Device;

/// The fastest device this build can use: CUDA, then Metal, then CPU.
///
/// GPU support is compile-time gated (`cuda` / `metal` features). A build
/// without them always returns the CPU, even on a machine with a GPU.
///
/// # Errors
/// Propagates the backend's initialisation error when a GPU is compiled in
/// but cannot be opened.
pub fn best_device() -> Result<Device, candle_core::Error> {
    #[cfg(feature = "cuda")]
    {
        match Device::new_cuda(0) {
            Ok(d) => return Ok(d),
            Err(e) => {
                // A machine with the feature compiled in but no usable GPU is
                // common (laptop without NVIDIA hardware); fall back rather than fail.
                eprintln!("CUDA unavailable ({e}); using CPU");
            }
        }
    }
    #[cfg(feature = "metal")]
    {
        match Device::new_metal(0) {
            Ok(d) => return Ok(d),
            Err(e) => eprintln!("Metal unavailable ({e}); using CPU"),
        }
    }
    Ok(Device::Cpu)
}

/// Human-readable description of a device, for logs and the UI.
pub fn device_description(device: &Device) -> String {
    match device {
        Device::Cpu => "CPU".to_string(),
        Device::Cuda(_) => "CUDA GPU".to_string(),
        Device::Metal(_) => "Metal GPU".to_string(),
    }
}
