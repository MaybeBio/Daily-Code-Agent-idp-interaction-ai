use candle_core::{Module, Tensor};
use candle_nn::{linear, Linear, VarBuilder};

/// A fully connected stack with GELU activations and a linear output.
///
/// Used by every head so the architecture is described in one place. No
/// dropout: these heads are small and are regularised by early stopping on a
/// held-out split instead, which keeps training and inference identical.
#[derive(Debug)]
pub struct Mlp {
    layers: Vec<Linear>,
}

impl Mlp {
    /// Build an MLP with the given layer widths, e.g. `[681, 256, 64, 1]`.
    ///
    /// # Errors
    /// Propagates tensor errors when a weight is missing or mis-shaped.
    pub fn new(widths: &[usize], vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        let layers = widths
            .windows(2)
            .enumerate()
            .map(|(i, w)| linear(w[0], w[1], vb.pp(format!("layer{i}"))))
            .collect::<candle_core::Result<Vec<_>>>()?;
        Ok(Self { layers })
    }
}

impl Module for Mlp {
    fn forward(&self, input: &Tensor) -> candle_core::Result<Tensor> {
        let mut x = input.clone();
        for (i, layer) in self.layers.iter().enumerate() {
            x = layer.forward(&x)?;
            // Activation on every layer but the output.
            if i + 1 < self.layers.len() {
                x = x.gelu_erf()?;
            }
        }
        Ok(x)
    }
}
