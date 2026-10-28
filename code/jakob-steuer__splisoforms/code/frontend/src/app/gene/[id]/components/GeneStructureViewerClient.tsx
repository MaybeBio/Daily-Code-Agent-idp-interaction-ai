"use client";

import React, { useEffect, useRef, useState } from "react";
import { RefreshCw, Box, AlertTriangle, Home } from "lucide-react";
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

const BG_DARK  = Color(0x0f172a);
const BG_LIGHT = Color(0xfafafa);

// ─── Shared Color Palette for Overlay Distinctions ───────────────────────────
const PALETTE = [
  0x3b82f6, // Blue
  0xef4444, // Red
  0x10b981, // Emerald
  0xf59e0b, // Amber
  0x8b5cf6, // Violet
  0xec4899, // Pink
  0x06b6d4, // Cyan
];

async function updateColor(plugin: PluginUIContext, struct: any, colorValue: number, alpha: number = 1.0) {
  if (!struct) return;
  for (const comp of struct.components) {
     for (const r of comp.representations) {
        await plugin.build().to(r.cell).update((old: any) => ({
           ...old,
           type: { ...old.type, params: { ...old.type?.params, alpha, ignoreLight: false, transparent: alpha < 1.0 } },
           colorTheme: { name: "uniform", params: { value: colorValue } }
        })).commit();
     }
  }
}

interface SelectedIsoform {
  id: number;
  transcript_id: string;
  /** Readable display name (MET-201 / MET-206-NOV-a6d00a); null for unmatched novel rows. */
  transcript_name?: string | null;
  sequence?: string;
  cifUrl?: string;
}

interface Props {
  canonicalCifUrl?: string;
  canonicalSequence?: string;
  selectedIsoforms: SelectedIsoform[];
  hoveredResidues?: number[];
  selectedResidues?: number[];
}

export default function GeneStructureViewerClient({
  canonicalCifUrl,
  canonicalSequence,
  selectedIsoforms,
  hoveredResidues = [],
  selectedResidues = []
}: Props) {
  const parentRef = useRef<HTMLDivElement>(null);
  const pluginRef = useRef<PluginUIContext | null>(null);
  const initRef   = useRef(false);

  const { theme } = useTheme();

  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [rmsds, setRmsds] = useState<Record<string, number>>({});

  // ── Sync canvas background with theme ────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current;
    if (!plugin?.canvas3d) return;
    plugin.canvas3d.setProps({ renderer: {
      ...plugin.canvas3d.props.renderer,
      backgroundColor: theme === "dark" ? BG_DARK : BG_LIGHT,
      selectColor: Color(0xFF1493), highlightColor: Color(0x39FF14),
    }, marking: {
      ...plugin.canvas3d.props.marking,
      enabled: true,
      highlightEdgeColor: Color(0x39FF14),
      selectEdgeColor: Color(0xFF1493),
      highlightEdgeStrength: 2,
      selectEdgeStrength: 2,
      edgeScale: 2
    }});
  }, [theme]);

  // ── Highlight / Hover ─────────────────────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current; if (!plugin) return;
    setTimeout(() => {
      if (hoveredResidues.length) {
        // Highlight in canonical structure
        const currentStructs = plugin.managers.structure.hierarchy.current.structures;
        const canonStruct = currentStructs[currentStructs.length - 1]; // Because we load it last and it stays
        const data = canonStruct?.cell.obj?.data;
        if (!data) return;
        const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
          "residue-test": Q.core.set.has([Q.set(...hoveredResidues), Q.ammp("label_seq_id")]) }), data);
        const loci = StructureSelection.toLociWithSourceUnits(sel);
        plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
      } else { plugin.managers.interactivity.lociHighlights.clearHighlights(); }
    }, 50);
  }, [hoveredResidues]);

  // ── Selection ─────────────────────────────────────────────────────────────────
  useEffect(() => {
    const plugin = pluginRef.current; if (!plugin) return;
    setTimeout(() => {
      const currentStructs = plugin.managers.structure.hierarchy.current.structures;
      const canonStruct = currentStructs[currentStructs.length - 1]; // Because we load it last and it stays
      const data = canonStruct?.cell.obj?.data;
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

  // 1. Initialize Mol* and Load Canonical Reference
  useEffect(() => {
    let alive = true;
    (async () => {
       if (!parentRef.current || initRef.current) return;
       initRef.current = true;

       const spec = DefaultPluginUISpec();
       spec.behaviors.push(PluginSpec.Behavior(MAQualityAssessment));
       spec.layout = { initial: { isExpanded: false, showControls: false, regionState: { left: "hidden", top: "hidden", right: "hidden", bottom: "hidden" } } };
       spec.components = { remoteState: "none" };

       const container = document.createElement("div");
       container.style.cssText = "width:100%;height:100%;position:absolute;overflow:hidden;touch-action:none;";
       parentRef.current.appendChild(container);

       let plugin: PluginUIContext;
       try { plugin = await createPluginUI({ target: container, render: renderReact18, spec }); }
       catch (e) { console.error("Mol* init:", e); if(alive) setLoading(false); return; }
       if (!alive) { plugin.dispose(); return; }
       pluginRef.current = plugin;

       // Theme Canvas
       plugin.canvas3d?.setProps({ renderer: { ...plugin.canvas3d.props.renderer, backgroundColor: theme === "dark" ? BG_DARK : BG_LIGHT } });

       // Load Canonical Core (Always at Index 0)
       if (canonicalCifUrl) {
           try {
               const BASE = getApiUrl();
               const url  = canonicalCifUrl.startsWith("http") ? canonicalCifUrl : `${BASE}${canonicalCifUrl}`;
               const dl = await plugin.builders.data.download({ url, isBinary: false, label: "Canonical" });
               const traj = await plugin.builders.structure.parseTrajectory(dl, "mmcif");
               await plugin.builders.structure.hierarchy.applyPreset(traj, "default", {
                   representationPreset: "polymer-and-ligand"
               });

               const structs = plugin.managers.structure.hierarchy.current.structures;
               if (structs.length > 0 && structs[0]?.cell.obj?.data) {
                  await updateColor(plugin, structs[0], 0x8899bb, 0.70); // Canonical — slate-silver, semi-transparent
               }
           } catch(e: any) {
               console.error("Canonical Load Failed:", e);
               if(alive) setError("Failed to instantiate canonical reference model.");
           }
       }
       if (alive) setLoading(false);
    })();

    return () => {
       alive = false;
       pluginRef.current?.dispose();
       pluginRef.current = null;
       initRef.current = false;
       if (parentRef.current) parentRef.current.innerHTML = "";
    };
  }, [canonicalCifUrl]);

  // 2. Synchronize Selected Isoforms Overlays
  useEffect(() => {
    const plugin = pluginRef.current;
    if (!plugin || loading || !canonicalCifUrl || !canonicalSequence) return;

    (async () => {
        // Remove ALL structures EXCEPT the canonical (which is structs[0])
        const structs = plugin.managers.structure.hierarchy.current.structures;
        if (structs.length > 1) {
            await plugin.managers.structure.hierarchy.remove(structs.slice(1));
        }

        const newRmsds: Record<string, number> = {};

        // Overlay sequentially
        for (let i = 0; i < selectedIsoforms.length; i++) {
           const iso = selectedIsoforms[i];
           if (!iso.cifUrl) continue;

           try {
               const BASE = getApiUrl();
               const url  = iso.cifUrl.startsWith("http") ? iso.cifUrl : `${BASE}${iso.cifUrl}`;
               const dl = await plugin.builders.data.download({ url, isBinary: false, label: iso.transcript_id });
               const traj = await plugin.builders.structure.parseTrajectory(dl, "mmcif");
               await plugin.builders.structure.hierarchy.applyPreset(traj, "default", {
                   representationPreset: "polymer-and-ligand"
               });

               // Allow tick to register structure
               await new Promise(r => setTimeout(r, 60));

               const currentStructs = plugin.managers.structure.hierarchy.current.structures;
               const thisIsoStruct = currentStructs[currentStructs.length - 1]; // Latest appended
               const canonStruct = currentStructs[0];

               // Apply Color
               const colorId = PALETTE[i % PALETTE.length];
               await updateColor(plugin, thisIsoStruct, colorId, 1.0);

               // Align and Superpose Geometry onto Canonical
               if (thisIsoStruct?.cell.obj?.data && canonStruct?.cell.obj?.data && iso.sequence) {
                   // Use Needleman-Wunsch alignment: only superpose on exact-match residue pairs.
                   // loci1 = canonical (reference, stays fixed), loci2 = isoform (mobile, gets bTransform).
                   const { idxCanon, idxQuery } = getAlignedIndices(canonicalSequence, iso.sequence);
                   const caQ = (idx: number[]) => (Q: any) => Q.struct.generator.atomGroups({
                      "atom-test":    Q.core.rel.eq([Q.ammp("label_atom_id"), "CA"]),
                      "residue-test": idx.length > 0 ? Q.core.set.has([Q.set(...idx), Q.ammp("label_seq_id")]) : Q.core.logic.and([]),
                   });
                   const loci1 = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(caQ(idxCanon), canonStruct.cell.obj.data));
                   const loci2 = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(caQ(idxQuery), thisIsoStruct.cell.obj.data));

                   if (idxQuery.length > 0) {
                       const { alignAndSuperpose } = await import("molstar/lib/mol-model/structure/structure/util/superposition");
                       const res = alignAndSuperpose([loci1, loci2]);
                       if (res && res.length && isFinite(res[0].rmsd)) {
                           await plugin.build().to(thisIsoStruct.cell).insert(StateTransforms.Model.TransformStructureConformation, {
                               transform: { name: "matrix", params: { data: res[0].bTransform } },
                           }).commit();
                           newRmsds[iso.transcript_id] = Math.round(res[0].rmsd * 100) / 100;
                       }
                   }
               }
           } catch(e) {
               console.warn(`Failed rendering Isoform ${iso.transcript_id}`, e);
           }
        }
        setRmsds(newRmsds);
    })();
  }, [selectedIsoforms, loading, canonicalCifUrl, canonicalSequence]);

  const resetCamera = () => {
    pluginRef.current?.canvas3d?.requestCameraReset();
  };

  return (
    <div className="relative w-full h-[500px] flex rounded-xl border border-primary/10 overflow-hidden" style={{ background: "transparent" }}>

      {/* Mol* UI overrides — hide the black settings box, keep the right toolbar */}
      <style>{`
        /* Hide layout region toggles and the top-left black settings panel */
        .msp-layout-top-region-toggle,
        .msp-layout-bottom-region-toggle,
        .msp-layout-left-region-toggle,
        .msp-layout-right-region-toggle,
        .msp-layout-hide-button,
        .msp-layout-static-toggle,
        .msp-sequence-viewer-toggle,
        .msp-state-snapshot-viewport-controls,
        .msp-animation-viewport-controls,
        .msp-viewport-controls-top-left,
        .msp-viewport-controls-bottom-left,
        .msp-viewport-controls-bottom-right {
          display: none !important;
          visibility: hidden !important;
          pointer-events: none !important;
        }
        /* Canvas background */
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

      <div ref={parentRef} className="absolute inset-0 z-10" />

      {/* Control Overlay — vertical card panel, theme-aware */}
      <div className="absolute top-4 left-4 z-20 flex flex-col gap-2 pointer-events-none">
        <div className="flex items-center gap-1 bg-background/95 px-3 py-1.5 rounded-lg border border-border/60 backdrop-blur-sm shadow-sm pointer-events-auto">
          <h3 className="font-bold text-foreground text-sm flex-1">Structure Viewer</h3>
          <button
            onClick={resetCamera}
            title="Reset camera"
            className="p-1 rounded hover:bg-primary/20 text-muted-foreground hover:text-primary transition-colors"
          >
            <Home className="w-3.5 h-3.5" />
          </button>
        </div>

        {loading && (
           <div className="flex items-center gap-2 bg-background/95 border border-border/60 text-foreground px-3 py-1.5 rounded-lg text-xs backdrop-blur-sm">
             <RefreshCw className="w-3.5 h-3.5 animate-spin text-primary" />
             Loading Canvas...
           </div>
        )}

        {error && (
           <div className="flex flex-col gap-1 bg-destructive/10 border border-destructive/30 text-destructive px-3 py-2 rounded-lg text-xs backdrop-blur-md max-w-sm pointer-events-auto">
             <span className="font-semibold flex items-center gap-1"><AlertTriangle className="w-3.5 h-3.5" /> Render Error</span>
             {error}
           </div>
        )}

        {/* Canonical legend entry */}
        {!loading && canonicalCifUrl && (
            <div className="flex items-center gap-2 bg-background/95 border border-border/60 px-3 py-1.5 rounded-lg backdrop-blur-sm w-full pointer-events-auto shadow-sm">
               <div className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: "#8899bb" }} />
               <span className="text-xs font-mono font-medium text-muted-foreground">Canonical Core</span>
            </div>
        )}

        {/* Isoform legend entries */}
        {!loading && selectedIsoforms.map((iso, i) => {
           const colorHex = `#${PALETTE[i % PALETTE.length].toString(16).padStart(6, '0')}`;
           const rmsd = rmsds[iso.transcript_id];
           return (
               <div key={iso.id} className="flex items-center justify-between gap-4 bg-background/95 border border-border/60 px-3 py-1.5 rounded-lg backdrop-blur-sm min-w-[200px] pointer-events-auto shadow-sm">
                   <div className="flex items-center gap-2">
                       <div className="w-2.5 h-2.5 rounded-full shadow-lg shrink-0" style={{ backgroundColor: colorHex, boxShadow: `0 0 6px ${colorHex}90` }} />
                       <span className="text-xs font-mono font-medium text-foreground">{iso.transcript_name || iso.transcript_id}</span>
                   </div>
                   {rmsd !== undefined && (
                       <span className="text-[10px] text-muted-foreground font-mono">RMSD: {rmsd}Å</span>
                   )}
               </div>
           );
        })}

        {!loading && selectedIsoforms.length === 0 && (
           <div className="bg-background/95 border border-border/60 px-3 py-1.5 rounded-lg backdrop-blur-sm text-xs text-muted-foreground italic pointer-events-auto">
               Select isoforms from the table to overlay structure.
           </div>
        )}
      </div>
    </div>
  );
}
