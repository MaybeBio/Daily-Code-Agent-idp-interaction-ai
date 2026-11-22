//! Protein language models and the prediction heads trained on top of them.
//!
//! The only model implemented so far is [`esm2`], a port of Meta's ESM-2
//! (Lin et al., 2023) to candle. The port is checked against reference
//! activations captured from the reference `PyTorch` implementation; see
//! `tests/esm2_reference.rs`.
//!
//! Weights are not bundled. [`weights`] downloads them from the Hugging Face
//! Hub on first use and caches them under the user's cache directory.

#![cfg_attr(test, allow(clippy::unwrap_used))]

mod device;
mod error;
pub mod esm2;
pub mod heads;
pub mod weights;

pub use device::{best_device, device_description};
pub use error::MlError;
