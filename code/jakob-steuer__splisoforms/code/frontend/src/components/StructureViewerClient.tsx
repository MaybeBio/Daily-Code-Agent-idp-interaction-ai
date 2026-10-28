"use client";

import React, { useEffect, useRef, useState, useCallback } from "react";
import { Download, Crosshair, RefreshCw, EyeOff, Box, ChevronDown, Home, Camera, Maximize } from "lucide-react";
import { createPluginUI } from "molstar/lib/mol-plugin-ui";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { DefaultPluginUISpec } from "molstar/lib/mol-plugin-ui/spec";
import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
import { Script } from "molstar/lib/mol-script/script";
import { StructureSelection } from "molstar/lib/mol-model/structure";
import { PluginSpec } from "molstar/lib/mol-plugin/spec";
import { MAQualityAssessment } from "molstar/lib/extensions/model-archive/quality-assessment/behavior";
import { Color } from "molstar/lib/mol-util/color";
import { StateTransforms } from "molstar/lib/mol-plugin-state/transforms";
import { getAlignedIndices } from "@/utils/alignment";
import { useTheme } from "@/context/ThemeContext";
import { getApiUrl } from "@/lib/api";
import "molstar/lib/mol-plugin-ui/skin/dark.scss";

const BG_DARK  = Color(0x0f172a); // slate-950
const BG_LIGHT = Color(0xfafafa); // near-white

// ─── Representation Modes ─────────────────────────────────────────────────────
type RepMode = "cartoon" | "backbone" | "surface" | "spacefill" | "ball-and-stick";
const REP_MODES: { id: RepMode; label: string; type: string }[] = [
  { id: "cartoon",   label: "Cartoon",   type: "cartoon" },
  { id: "backbone",  label: "Backbone",  type: "backbone" },
  { id: "surface",   label: "Surface",   type: "molecular-surface" },
  { id: "spacefill", label: "Spacefill", type: "spacefill" },
  { id: "ball-and-stick", label: "Ball & Stick", type: "ball-and-stick" },
];


async function updateColor(plugin: any, struct: any, colorValue: number = 0x8899bb, alpha: number = 0.88) {
  if (!struct) return;

  for (const comp of struct.components) {
    // Keep any co-folded ligand (drug) solid — only ghost the protein.
    if (isLigandComp(comp)) continue;
    for (const r of comp.representations) {
      await plugin.build().to(r.cell).update((old: any) => ({
         ...old,
         type: {
            ...old.type,
            params: {
               ...old.type?.params,
               alpha,
               ignoreLight: true,
               transparent: alpha < 1.0
            }
         },
         colorTheme: {
            name: "uniform",
            params: { value: colorValue }
         }
      })).commit();
    }
  }
}

// ─── Props ────────────────────────────────────────────────────────────────────
interface Props {
  cifUrl: string;
  canonicalCifUrl?: string;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  junctionRange?: [number, number];
  sequence?: string;
  canonicalSequence?: string;
  onPlddtCalculated?: (plddt: number) => void;
  // ── Drug-complex opt-in behaviours (default off; isoform pages unaffected) ──
  ligandBallAndStick?: boolean;
  pocketResidues?: number[];
  /** Residues the "Show binding pocket" toggle highlights on the PRIMARY structure.
   *  When omitted, the toggle falls back to pocketResidues. On the isoform viewer the
   *  primary is the isoform, so pass the pocket in ISOFORM register (retained residues)
   *  here — pocketResidues stays in canonical register for focusPocketOnCanonical. */
  pocketResiduesPrimary?: number[];
  enablePocketToggle?: boolean;
  focusPocketOnCanonical?: boolean;   // on "show canonical": ball-and-stick + H-bonds on the pocket
  canonicalLabel?: string;
  /** Legend label for the structure being displayed. Defaults to "Novel Isoform", which is
   *  wrong on a canonical isoform's own page — pass "Canonical Isoform" there. */
  primaryLabel?: string;
  /** Show ONLY the co-folded ligand: hide the protein and frame the camera on the drug.
   *  Used by the drug summary panel, where the molecule itself is the subject. The
   *  coordinates are the real bound conformer from our own CIF — not a generic one. */
  ligandOnly?: boolean;
}

// Is a structure component the (co-folded) ligand rather than the protein polymer?
const isPolymerComp = (c: any) =>
  c?.key === "polymer" || /polymer/i.test(c?.cell?.obj?.label || "");
const isLigandComp = (c: any) =>
  c?.key === "ligand" || /ligand/i.test(c?.cell?.obj?.label || "");

export default function StructureViewerClient({
  cifUrl, canonicalCifUrl,
  hoveredResidues = [], selectedResidues = [],
  junctionRange, sequence, canonicalSequence,
  onPlddtCalculated,
  ligandBallAndStick = false, pocketResidues = [],
  pocketResiduesPrimary,
  enablePocketToggle = false, focusPocketOnCanonical = false,
  canonicalLabel,
  primaryLabel = "Novel Isoform",
  ligandOnly = false,
}: Props) {
  if (typeof window === "undefined") return null;

  const parentRef   = useRef<HTMLDivElement>(null);
  const pluginRef   = useRef<PluginUIContext | null>(null);
  const initRef     = useRef(false);
  // Track whether canonical is currently loaded in the state tree
  const canonLoaded = useRef(false);

  const { theme } = useTheme();

  const [error,       setError]       = useState<string | null>(null);
  const [loading,     setLoading]     = useState(true);
  const [showCanon,   setShowCanon]   = useState(false);
  const [canonBusy,   setCanonBusy]   = useState(false);  // loading / removing
  const [rmsd,        setRmsd]        = useState<number | null>(null);
  const [alignInfo,   setAlignInfo]   = useState<{ n: number; pocket: boolean } | null>(null);
  const [repMode,     setRepMode]     = useState<RepMode>("cartoon");
  const [showRepMenu, setShowRepMenu] = useState(false);
  const [showPocket,  setShowPocket]  = useState(false);

  // ── Sync canvas background with theme ────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current;
    if (!plugin?.canvas3d) return;
    plugin.canvas3d.setProps({ renderer: {
      ...plugin.canvas3d.props.renderer,
      backgroundColor: theme === "dark" ? BG_DARK : BG_LIGHT,
    }});
  }, [theme]);

  // ── Plugin + novel structure init ────────────────────────────────────────────
  useEffect(() => {
    let alive = true;
    (async () => {
      if (!alive || !parentRef.current || initRef.current) return;
      initRef.current = true;

      const spec = DefaultPluginUISpec();
      spec.behaviors.push(PluginSpec.Behavior(MAQualityAssessment));
      spec.layout = { initial: {
        isExpanded: false, showControls: false,
        regionState: { left: "hidden", top: "hidden", right: "hidden", bottom: "hidden" },
      }};
      spec.components = { 
        remoteState: "none",
        viewport: {
          controls: () => null // Hide the faint right-side toolbar
        }
      };

      const container = document.createElement("div");
      container.style.cssText = "width:100%;height:100%;position:absolute;overflow:hidden;touch-action:none;";
      parentRef.current!.appendChild(container);

      let plugin: PluginUIContext;
      try { plugin = await createPluginUI({ target: container, render: renderReact18, spec }); }
      catch (e) { console.error("Mol* init:", e); if (alive) setLoading(false); return; }
      if (!alive) { plugin.dispose(); return; }
      pluginRef.current = plugin;

      // Canvas background (theme-aware) + custom selection/highlight colours
      plugin.canvas3d?.setProps({ renderer: {
        ...plugin.canvas3d.props.renderer,
        backgroundColor: theme === "dark" ? BG_DARK : BG_LIGHT,
        selectColor: Color(0xFF1493), highlightColor: Color(0x39FF14),
      }});

      // Load novel structure with pLDDT
      try {
        const BASE = getApiUrl();
        const url  = cifUrl.startsWith("http") ? cifUrl : `${BASE}${cifUrl}`;
        const dl   = await plugin.builders.data.download({ url, isBinary: false });
        const traj = await plugin.builders.structure.parseTrajectory(dl, "mmcif");
        await plugin.builders.structure.hierarchy.applyPreset(traj, "default", {
          representationPreset: "polymer-and-ligand",
          representationPresetParams: { theme: { globalName: "plddt-confidence" } },
        });

        // Ligand-only: drop the protein components entirely (hiding them still lets the
        // camera auto-frame include them, so the drug ends up a speck) and frame the drug.
        if (ligandOnly) {
          const struct = plugin.managers.structure.hierarchy.current.structures[0];
          const polymers = (struct?.components ?? []).filter(isPolymerComp);
          if (polymers.length) {
            await plugin.managers.structure.hierarchy.remove(polymers);
          }

          // Colour the drug by ELEMENT (CPK), not pLDDT. The preset applies
          // plddt-confidence globally, which paints the whole molecule one flat confidence
          // colour — useless for reading a small molecule, where the point is which atom is
          // N vs O vs S. AF3 CIFs carry heavy atoms only (no hydrogens at all), so there is
          // nothing to add there without generating coordinates via a chemistry toolkit.
          const ligands = (plugin.managers.structure.hierarchy.current.structures[0]?.components ?? [])
            .filter(isLigandComp);
          const b = plugin.build();
          for (const c of ligands) {
            for (const r of c.representations) {
              b.to(r.cell).update(
                StateTransforms.Representation.StructureRepresentation3D,
                (old: any) => ({
                  ...old,
                  colorTheme: { name: "element-symbol", params: {} },
                }),
              );
            }
          }
          await b.commit();

          // Drop the XYZ axes widget — orientation gizmo clutter on a thumbnail.
          plugin.canvas3d?.setProps({
            camera: { helper: { axes: { name: "off", params: {} } } },
          });

          plugin.managers.camera.reset();
        }

        // Report avg pLDDT
        const mdl = plugin.managers.structure.hierarchy.current.models[0]?.cell.obj?.data;
        if (mdl && onPlddtCalculated) {
          const bf = Array.from(mdl.atomicConformation.B_iso_or_equiv.toArray() || []);
          if (bf.length) onPlddtCalculated(bf.reduce((s, v) => s + v, 0) / bf.length);
        }
      } catch (e: any) { 
        console.warn("Novel CIF load failed:", e); 
        if (alive) setError(e.message || "Failed to load structure model.");
      }

      if (alive) setLoading(false);
    })();

    return () => {
      alive = false;
      pluginRef.current?.dispose();
      pluginRef.current = null;
      initRef.current = false;
      canonLoaded.current = false;
      if (parentRef.current) parentRef.current.innerHTML = "";
    };
  }, [cifUrl, ligandOnly]);

  // ── Canonical toggle (load-on-click, remove via hierarchy manager) ────────────
  const handleToggleCanonical = useCallback(async () => {
    const plugin = pluginRef.current;
    if (!plugin || canonBusy) return;
    setCanonBusy(true);

    try {
      if (canonLoaded.current) {
        // ── REMOVE ───────────────────────────────────────────────────────────
        // Use the official hierarchy manager remove() — removes all non-first structures
        const structs = plugin.managers.structure.hierarchy.current.structures;
        if (structs.length > 1) {
          // remove() accepts StructureRef[], toggling their visibility in the hierarchy
          await plugin.managers.structure.hierarchy.remove(structs.slice(1));
        }
        // Clear the pocket focus (ball-and-stick + H-bonds) view.
        if (focusPocketOnCanonical) {
          plugin.managers.structure.focus.clear();
        }
        canonLoaded.current = false;
        setShowCanon(false);
        setRmsd(null);
        setAlignInfo(null);

      } else {
        // ── LOAD ─────────────────────────────────────────────────────────────
        if (!canonicalCifUrl) return;
        const BASE = getApiUrl();
        const url  = canonicalCifUrl.startsWith("http") ? canonicalCifUrl : `${BASE}${canonicalCifUrl}`;

        console.log("[Canonical] Loading:", url);

        const dl   = await plugin.builders.data.download({ 
          url, 
          isBinary: false, 
          label: `Canonical Reference (${Date.now()})` 
        });
        const traj = await plugin.builders.structure.parseTrajectory(dl, "mmcif");

        // Use default preset (no pLDDT)
        await plugin.builders.structure.hierarchy.applyPreset(traj, "default", {
          representationPreset: "polymer-and-ligand",
        });

        canonLoaded.current = true;

        // Give the state tree a tick to settle, then recolour white
        await new Promise(r => setTimeout(r, 80));

        const structs = plugin.managers.structure.hierarchy.current.structures;
        const cIdx = structs.length - 1;   // canonical is always the last-loaded structure
        console.log("[Canonical] Structure count:", structs.length, "using index", cIdx);

        if (cIdx > 0) {
          setShowCanon(true);

          // ── Apply styling ────────────────────────────────────────────────
          // Ghost the canonical protein (silver); the novel isoform stays solid
          // (pLDDT). On drug pages the co-folded ligand is kept solid by
          // updateColor, so the drug stays prominent over the ghosted protein.
          const canonStruct = structs[cIdx];
          const novelStruct = structs[0];
          if (canonStruct?.cell.obj?.data) await updateColor(plugin, canonStruct);
          if (novelStruct?.cell.obj?.data && canonStruct?.cell.obj?.data) {
            // Splice-aware alignment → exact-match shared-exon Cα pairs (global anchors).
            const { idxCanon, idxQuery } = getAlignedIndices(canonicalSequence || '', sequence || '');

            // Prefer superposing on the RETAINED BINDING-POCKET residues so the
            // pocket overlays cleanly. A global fit averages over the whole fold
            // and leaves the pocket misaligned when domains shift between the two
            // independently-predicted AF3 models. Fall back to the global set when
            // too little of the pocket survives (≤25%) or too few anchors remain.
            let supCanon = idxCanon, supQuery = idxQuery, pocketAligned = false;
            if (focusPocketOnCanonical && pocketResidues.length) {
              const canonToIso = new Map<number, number>();
              for (let k = 0; k < idxCanon.length; k++) canonToIso.set(idxCanon[k], idxQuery[k]);
              const pc: number[] = [], pq: number[] = [];
              for (const p of pocketResidues) {
                const iso = canonToIso.get(p);
                if (iso !== undefined) { pc.push(p); pq.push(iso); }
              }
              if (pc.length >= 3 && pc.length > 0.25 * pocketResidues.length) {
                supCanon = pc; supQuery = pq; pocketAligned = true;
              }
            }

            const caQ = (idx: number[]) => (Q: any) => Q.struct.generator.atomGroups({
              "atom-test":    Q.core.rel.eq([Q.ammp("label_atom_id"), "CA"]),
              "residue-test": idx.length > 0
                ? Q.core.set.has([Q.set(...idx), Q.ammp("label_seq_id")])
                : Q.core.logic.and([]),
            });
            const loci1 = StructureSelection.toLociWithSourceUnits(
              Script.getStructureSelection(caQ(supQuery), novelStruct.cell.obj.data));
            const loci2 = StructureSelection.toLociWithSourceUnits(
              Script.getStructureSelection(caQ(supCanon), canonStruct.cell.obj.data));
            try {
              console.log("[Superposition]", pocketAligned ? "pocket-local" : "global",
                          { n: supCanon.length });
              const { alignAndSuperpose } = await import("molstar/lib/mol-model/structure/structure/util/superposition");
              const res = alignAndSuperpose([loci1, loci2]);
              if (res?.length) {
                await plugin.build()
                  .to(canonStruct.cell)
                  .insert(StateTransforms.Model.TransformStructureConformation, {
                    transform: { name: "matrix", params: { data: res[0].bTransform } },
                  }).commit();
                if (isFinite(res[0].rmsd)) {
                  setRmsd(Math.round(res[0].rmsd * 100) / 100);
                  setAlignInfo({ n: supCanon.length, pocket: pocketAligned });
                }
              }
            } catch (e) { console.warn("Superposition:", e); }
          }

          // Show the drug pocket as ball-and-stick with interaction/H-bond dashes
          // — the same view Mol* gives on a residue click. Focus residues sit on
          // the canonical (which carries the ligand), rendered over the ghost.
          if (focusPocketOnCanonical && pocketResidues.length && canonStruct?.cell.obj?.data) {
            try {
              const cdata = canonStruct.cell.obj.data;
              const psel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
                "residue-test": Q.core.set.has([Q.set(...pocketResidues), Q.ammp("label_seq_id")]),
              }), cdata);
              const ploci = StructureSelection.toLociWithSourceUnits(psel);
              plugin.managers.structure.focus.setFromLoci(ploci);  // ball-and-stick + H-bonds
              plugin.managers.camera.focusLoci(ploci);             // zoom the camera to the pocket
            } catch (e) { console.warn("pocket focus:", e); }
          }
        }
      }
    } catch (e) {
      console.error("[Canonical] toggle failed:", e);
      canonLoaded.current = false;
      setShowCanon(false);
    } finally {
      setCanonBusy(false);
    }
  }, [canonBusy, canonicalCifUrl, sequence, canonicalSequence, focusPocketOnCanonical, pocketResidues]);

  // ── Representation mode ───────────────────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current; if (!plugin || loading) return;
    const type = REP_MODES.find(m => m.id === repMode)?.type ?? "cartoon";
    const structs = plugin.managers.structure.hierarchy.current.structures;
    const novel   = structs[0]; if (!novel) return;
    const b = plugin.build();
    for (const c of novel.components) {
      // Keep the co-folded drug in ball&stick; the switcher only restyles the protein.
      if (ligandBallAndStick && isLigandComp(c)) continue;
      for (const r of c.representations)
        b.to(r.cell).update(StateTransforms.Representation.StructureRepresentation3D,
          (old: any) => ({ ...old, type: { name: type, params: old.type?.params ?? {} } }));
    }
    b.commit().catch(console.warn);
  }, [repMode, loading, ligandBallAndStick]);

  // ── Hover ─────────────────────────────────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current; if (!plugin) return;
    setTimeout(() => {
      const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!data) return;
      if (hoveredResidues.length) {
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
          "residue-test": Q.core.set.has([Q.set(...hoveredResidues), Q.ammp("label_seq_id")]) }), data);
        plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: StructureSelection.toLociWithSourceUnits(sel) });
      } else { plugin.managers.interactivity.lociHighlights.clearHighlights(); }
    }, 50);
  }, [hoveredResidues]);

  // ── Selection ─────────────────────────────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current; if (!plugin) return;
    setTimeout(() => {
      const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!data) return;
      if (selectedResidues.length > 0) {
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
          "residue-test": Q.core.set.has([Q.set(...selectedResidues), Q.ammp("label_seq_id")]) }), data);
        const loci = StructureSelection.toLociWithSourceUnits(sel);
        plugin.managers.camera.focusLoci(loci);
        plugin.managers.interactivity.lociSelects.selectOnly({ loci });
      } else {
        plugin.managers.interactivity.lociSelects.deselectAll();
      }
    }, 50);
  }, [selectedResidues]);

  // ── Binding-pocket toggle (drug pages) ─────────────────────────────────────────
  // The toggle highlights residues on the PRIMARY structure, so it must use the
  // pocket in that structure's register: pocketResiduesPrimary (isoform register on
  // the isoform viewer) when supplied, else pocketResidues (canonical register).
  const togglePocket = (pocketResiduesPrimary ?? pocketResidues);
  useEffect(() => {
    if (!enablePocketToggle) return;
    const plugin = pluginRef.current; if (!plugin || loading) return;
    setTimeout(() => {
      const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
      if (!data) return;
      if (showPocket && togglePocket.length > 0) {
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
          "residue-test": Q.core.set.has([Q.set(...togglePocket), Q.ammp("label_seq_id")]) }), data);
        const loci = StructureSelection.toLociWithSourceUnits(sel);
        plugin.managers.camera.focusLoci(loci);
        plugin.managers.interactivity.lociSelects.selectOnly({ loci });
      } else {
        plugin.managers.interactivity.lociSelects.deselectAll();
      }
    }, 50);
  }, [showPocket, togglePocket, enablePocketToggle, loading]);

  // ── Spin / Focus ──────────────────────────────────────────────────────────────
  const handleSpin = () => {
    const c3d = pluginRef.current?.canvas3d; if (!c3d) return;
    const on = c3d.props.trackball.animate?.name === "spin";
    c3d.setProps({ trackball: { animate: on ? { name: "off", params: {} } : { name: "spin", params: { speed: 1, axis: [0, 1, 0] } } } });
  };
  const handleFocus = () => {
    const plugin = pluginRef.current; if (!plugin || !junctionRange) return;
    const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
    if (!data) return;
    const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
      "residue-test": Q.core.rel.inRange([Q.ammp("label_seq_id"), junctionRange[0], junctionRange[1]]) }), data);
    plugin.managers.camera.focusLoci(StructureSelection.toLociWithSourceUnits(sel));
  };

  const repLabel = REP_MODES.find(m => m.id === repMode)?.label ?? "Cartoon";

  // ─────────────────────────────────────────────────────────────────────────────
  const ICON_BTN = "p-1.5 rounded hover:bg-muted hover:text-foreground transition text-muted-foreground";

  return (
    <div className="relative w-full h-full flex flex-col rounded-xl overflow-hidden border border-primary/10" style={{ background: "transparent" }}>

      {/* Loading overlay */}
      {loading && (
        <div className="absolute inset-0 z-50 flex items-center justify-center bg-background/70 backdrop-blur-sm">
           <div className="flex flex-col items-center gap-3">
             <RefreshCw className="w-8 h-8 text-primary animate-spin" />
             <span className="text-sm font-medium text-muted-foreground">Loading 3D Model...</span>
           </div>
        </div>
      )}

      {/* Error overlay */}
      {error && (
        <div className="absolute inset-0 z-50 flex items-center justify-center bg-background/90 backdrop-blur-md p-6">
           <div className="flex flex-col items-center gap-4 max-w-md text-center">
             <div className="w-12 h-12 rounded-full bg-destructive/10 flex items-center justify-center border border-destructive/20">
               <Box className="w-6 h-6 text-destructive" />
             </div>
             <div className="flex flex-col gap-1">
               <h4 className="text-lg font-semibold text-foreground">Model Load Failed</h4>
               <p className="text-sm text-muted-foreground break-words">{error}</p>
             </div>
             <button onClick={() => window.location.reload()} className="px-4 py-2 bg-muted hover:bg-muted/80 text-foreground rounded-md text-sm transition-colors border border-border">
               Retry Loading
             </button>
           </div>
        </div>
      )}

      {/* Mol* canvas */}
      <style>{`
        /* Hide Molstar default UI panels and toggles */
        .msp-layout-top-region-toggle,
        .msp-layout-bottom-region-toggle,
        .msp-layout-left-region-toggle,
        .msp-layout-right-region-toggle,
        .msp-layout-hide-button,
        .msp-layout-static-toggle,
        .msp-sequence-viewer-toggle,
        .msp-viewport-controls,
        .msp-state-snapshot-viewport-controls,
        .msp-animation-viewport-controls,
        .msp-viewport-controls-bottom-left,
        .msp-viewport-controls-bottom-right,
        .msp-viewport-controls-top-left {
          display: none !important;
          visibility: hidden !important;
          opacity: 0 !important;
          pointer-events: none !important;
        }
        
        /* Ensure the Molstar canvas container itself is visible */
        .msp-plugin-content {
          background: transparent !important;
        }
        /* Kill Mol*'s default canvas/plugin borders + any pseudo-element
           lines + descendant borders. Belt and braces — the previous
           targeted reset still left a faint 1px line at the corners on
           light-mode. */
        .msp-plugin,
        .msp-plugin *,
        .msp-plugin canvas,
        .msp-plugin-content,
        .msp-plugin-content *,
        .msp-viewport,
        .msp-viewport * {
          border: 0 !important;
          outline: 0 !important;
          box-shadow: none !important;
        }
        .msp-plugin::before, .msp-plugin::after,
        .msp-plugin-content::before, .msp-plugin-content::after,
        .msp-viewport::before, .msp-viewport::after {
          display: none !important;
        }
      `}</style>
      <div ref={parentRef} className="absolute inset-0" />

      {/* ── Left panel — vertical cards, mirrors gene viewer style ──
           In ligandOnly mode this whole stack (title bar, rep switcher, legend chip) is
           dropped: at thumbnail size it obscures the molecule it is meant to describe. ── */}
      <div className={`absolute top-4 left-4 z-20 flex-col gap-2 pointer-events-none ${ligandOnly ? "hidden" : "flex"}`}>

        {/* Title card + icon controls inline — z-10 keeps its stacking context above sibling cards so the dropdown renders in front */}
        <div className="relative z-10 bg-background/95 border border-border/60 px-3 py-1.5 rounded-lg backdrop-blur-sm shadow-sm pointer-events-auto flex items-center justify-between gap-3">
          <h3 className="font-bold text-sm text-foreground whitespace-nowrap">Structure Viewer</h3>
          <div className="flex items-center gap-0.5">
            {/* Representation dropdown */}
            <div className="relative">
              <button title={`Representation: ${repLabel}`} onClick={() => setShowRepMenu(v => !v)} className={ICON_BTN}>
                <Box className="w-3.5 h-3.5" />
              </button>
              {showRepMenu && (
                <div className="absolute left-0 top-full mt-1 bg-background border border-border rounded-lg shadow-xl z-50 min-w-[120px] overflow-hidden">
                  {REP_MODES.map(m => (
                    <button key={m.id} onClick={() => { setRepMode(m.id); setShowRepMenu(false); }}
                      className={`block w-full text-left px-3 py-1.5 text-[11px] hover:bg-muted transition ${repMode === m.id ? "text-primary font-semibold" : "text-foreground/70"}`}>
                      {m.label}
                    </button>
                  ))}
                </div>
              )}
            </div>
            <div className="w-px h-3.5 bg-border mx-0.5" />
            <button title="Spin" onClick={handleSpin} className={ICON_BTN}>
              <RefreshCw className="w-3.5 h-3.5" />
            </button>
            <button title="Reset view" onClick={() => pluginRef.current?.managers.camera.reset()} className={ICON_BTN}>
              <Home className="w-3.5 h-3.5" />
            </button>
            <button title="Take screenshot" onClick={() => pluginRef.current?.helpers.viewportScreenshot?.download('structure.png')} className={ICON_BTN}>
              <Camera className="w-3.5 h-3.5" />
            </button>
            <button title="Full screen" onClick={() => {
              if (!document.fullscreenElement) {
                parentRef.current?.requestFullscreen().catch(console.error);
              } else {
                document.exitFullscreen();
              }
            }} className={ICON_BTN}>
              <Maximize className="w-3.5 h-3.5" />
            </button>
            <button title="Download CIF" onClick={() => window.open(cifUrl, "_blank")} className={ICON_BTN}>
              <Download className="w-3.5 h-3.5" />
            </button>
          </div>
        </div>

        {/* Displayed-structure legend entry */}
        {!loading && (
          <div className="self-start flex items-center gap-2 bg-background/95 border border-border/60 px-3 py-1.5 rounded-lg backdrop-blur-sm shadow-sm">
            <div className="w-2.5 h-2.5 rounded-full shrink-0 bg-gradient-to-br from-[#0053D6] to-[#FF7D45]" />
            <span className="text-xs font-mono font-medium text-foreground">{primaryLabel}</span>
            <span className="text-[10px] text-muted-foreground">(pLDDT)</span>
          </div>
        )}

        {/* Binding-pocket toggle (drug pages) */}
        {!loading && enablePocketToggle && (
          <button
            onClick={() => setShowPocket(v => !v)}
            className={`self-start flex items-center gap-2 px-3 py-1.5 rounded-lg backdrop-blur-sm shadow-sm border pointer-events-auto transition-colors text-xs font-mono font-medium
              ${showPocket
                ? "bg-[#FF1493]/10 border-[#FF1493]/40 text-[#FF1493]"
                : "bg-background/95 border-border/60 text-foreground hover:bg-muted/60"}`}
            title="Highlight drug-contact residues (later: full binding-pocket analysis)"
          >
            {showPocket
              ? <div className="w-2.5 h-2.5 rounded-full shrink-0" style={{ backgroundColor: "#FF1493" }} />
              : <Crosshair className="w-3 h-3 shrink-0 opacity-60" />}
            <span>{showPocket ? `Binding pocket (${togglePocket.length})` : "Show binding pocket"}</span>
          </button>
        )}

        {/* Canonical toggle — styled like a legend entry when active */}
        {!loading && canonicalCifUrl && (
          <button
            onClick={handleToggleCanonical}
            disabled={canonBusy}
            className={`self-start flex items-center gap-2 px-3 py-1.5 rounded-lg backdrop-blur-sm shadow-sm border pointer-events-auto transition-colors text-xs font-mono font-medium
              ${showCanon
                ? "bg-primary/10 border-primary/30 text-primary"
                : "bg-background/95 border-border/60 text-foreground hover:bg-muted/60"}`}
          >
            {canonBusy
              ? <RefreshCw className="w-3 h-3 animate-spin shrink-0" />
              : showCanon
                ? <div className="w-2.5 h-2.5 rounded-full shrink-0" style={{ backgroundColor: "#8899bb" }} />
                : <EyeOff className="w-3 h-3 shrink-0 opacity-50" />
            }
            <span>{canonBusy ? "Loading…" : showCanon ? (canonicalLabel ? `${canonicalLabel} (silver)` : "Canonical (silver)") : (canonicalLabel ? `Show ${canonicalLabel}` : "Show Canonical")}</span>
            {showCanon && rmsd !== null && alignInfo && (
              <span className="text-[10px] font-normal ml-auto text-muted-foreground whitespace-nowrap">
                RMSD: {rmsd.toFixed(2)} Å · {alignInfo.n} {alignInfo.pocket ? "pocket Cα" : "Cα"}
              </span>
            )}
          </button>
        )}
      </div>

      {/* ── pLDDT legend (bottom-right) — theme-aware ──
           Hidden in ligandOnly mode: a small molecule has no meaningful per-residue
           confidence, and the legend would cover most of a thumbnail-sized viewport. ── */}
      <div className={`absolute bottom-3 right-3 z-10 bg-background/95 backdrop-blur-sm p-2.5 rounded-lg border border-border/60 shadow-sm pointer-events-none ${ligandOnly ? "hidden" : ""}`}>
        <div className="text-[10px] font-semibold text-muted-foreground mb-1.5 uppercase tracking-wider">pLDDT</div>
        <div className="flex flex-col gap-0.5 text-[10px]">
          {[["#0053D6",">90 Very High"],["#65CBF3","70–90 High"],["#FFDB13","50–70 Low"],["#FF7D45","<50 Very Low"]].map(([c,l]) => (
            <div key={l} className="flex items-center gap-1.5">
              <div className="w-2.5 h-2.5 rounded-full shrink-0" style={{ backgroundColor: c }} />
              <span className="text-muted-foreground">{l}</span>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
