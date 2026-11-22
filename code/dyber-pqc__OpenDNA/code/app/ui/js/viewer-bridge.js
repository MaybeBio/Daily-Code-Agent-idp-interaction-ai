// Thin wrapper around the Mol* viewer, called from Rust via wasm-bindgen.
// Keep this file small: logic belongs in Rust, this only adapts Mol*'s API.
"use strict";

window.OpenDNAViewer = (() => {
  let viewer = null;
  let host = null;
  let resizeObserver = null;

  // Match the canvas to the app theme (CSS token --surface-sunken).
  function applyBackground() {
    if (!viewer || !viewer.plugin.canvas3d) return;
    const css = getComputedStyle(document.documentElement).getPropertyValue("--surface-sunken").trim();
    const color = parseInt(css.replace("#", ""), 16);
    if (!Number.isNaN(color)) viewer.plugin.canvas3d.setProps({ renderer: { backgroundColor: color } });
  }
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", applyBackground);

  async function ensure(element) {
    if (viewer && host === element) return viewer;
    if (viewer) viewer.dispose();
    if (resizeObserver) resizeObserver.disconnect();
    host = element;
    viewer = await molstar.Viewer.create(element, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      collapseLeftPanel: true,
      viewportShowExpand: false,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false,
      viewportShowTrajectoryControls: false,
      volumeStreamingDisabled: true,
    });
    applyBackground();
    resizeObserver = new ResizeObserver(() => viewer && viewer.handleResize());
    resizeObserver.observe(element);
    return viewer;
  }

  return {
    async load(element, content, format) {
      const v = await ensure(element);
      await v.plugin.clear();
      v.handleResize();
      await v.loadStructureFromData(content, format === "mmcif" ? "mmcif" : "pdb");
      v.plugin.managers.camera.reset();
    },
    clear() {
      if (viewer) viewer.plugin.clear();
    },
    highlight(chain, residue) {
      if (!viewer) return;
      viewer.structureInteractivity({
        action: "highlight",
        elements: residue === null ? undefined : { auth_asym_id: chain, auth_seq_id: residue },
      });
    },
    focus(chain, residue) {
      if (!viewer) return;
      viewer.structureInteractivity({ action: "focus", elements: { auth_asym_id: chain, auth_seq_id: residue } });
    },
  };
})();
