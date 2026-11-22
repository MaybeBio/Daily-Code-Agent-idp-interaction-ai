//! Bindings to `js/viewer-bridge.js`, the only JavaScript in the app.

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = OpenDNAViewer, js_name = load, catch)]
    async fn js_load(host: &web_sys::HtmlElement, content: &str, format: &str) -> Result<JsValue, JsValue>;
    #[wasm_bindgen(js_namespace = OpenDNAViewer, js_name = clear)]
    pub fn clear();
    #[wasm_bindgen(js_namespace = OpenDNAViewer, js_name = highlight)]
    fn js_highlight(chain: &str, residue: JsValue);
    #[wasm_bindgen(js_namespace = OpenDNAViewer, js_name = focus)]
    pub fn focus(chain: &str, residue: i32);
}

/// Load a coordinate file into the viewer inside `host`.
pub async fn load(host: &web_sys::HtmlElement, content: &str, format: &str) -> Result<(), String> {
    js_load(host, content, format).await.map(|_| ()).map_err(|e| format!("viewer: {e:?}"))
}

/// Highlight one residue, or clear the highlight with `None`.
pub fn highlight(chain: &str, residue: Option<i32>) {
    js_highlight(chain, residue.map_or(JsValue::NULL, JsValue::from));
}
