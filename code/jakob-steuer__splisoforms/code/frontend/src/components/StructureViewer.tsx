"use client";

import dynamic from "next/dynamic";
import { RefreshCw } from "lucide-react";

const StructureViewerClient = dynamic(() => import("./StructureViewerClient"), {
  ssr: false,
  loading: () => (
    <div className="relative w-full h-full min-h-[400px] bg-slate-900/50 rounded-lg overflow-hidden border border-slate-800 flex items-center justify-center">
       <div className="flex flex-col items-center gap-3">
         <RefreshCw className="w-8 h-8 text-primary animate-spin" />
         <span className="text-sm font-medium text-muted-foreground">Loading Structure Viewer...</span>
       </div>
    </div>
  ),
});

interface ViewerProps {
  cifUrl: string;
  canonicalCifUrl?: string;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  junctionRange?: [number, number];
  sequence?: string;
  canonicalSequence?: string;
  onPlddtCalculated?: (plddt: number) => void;
  // ── Drug-complex opt-in behaviours (default off; isoform pages unaffected) ──
  ligandBallAndStick?: boolean;      // keep ligand ball&stick; rep switcher only changes the protein
  pocketResidues?: number[];         // residues for the "Show binding pockets" toggle
  pocketResiduesPrimary?: number[];  // pocket in the PRIMARY structure's register (isoform viewer)
  enablePocketToggle?: boolean;      // render the pocket toggle button
  focusPocketOnCanonical?: boolean;  // on "Show canonical": ball-and-stick + H-bonds on the pocket
  canonicalLabel?: string;           // override the canonical toggle label
  primaryLabel?: string;             // legend label for the displayed structure (default "Novel Isoform")
  ligandOnly?: boolean;              // show only the co-folded ligand (drug summary panel)
}

export default function StructureViewer(props: ViewerProps) {
  return <StructureViewerClient {...props} />;
}
