//! Hardware detection used to pick model sizes and compute backends.

#![cfg_attr(test, allow(clippy::unwrap_used))]
mod detect;

pub use detect::{
    detect_hardware, get_device_capability, ComputeBackend, DeviceCapability, GpuInfo, HardwareInfo,
    ModelTier, Precision,
};
