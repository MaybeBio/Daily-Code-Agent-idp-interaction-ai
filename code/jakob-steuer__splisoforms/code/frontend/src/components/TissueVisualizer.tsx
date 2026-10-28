"use client";

import React, { useMemo } from "react";
import { cn } from "@/lib/utils";

type OrganCoordinate = {
  x: number;
  y: number;
  labelY: number; // custom Y for label to compactly fit within 0..250 height
  label: string;
  side: "left" | "right";
};

// Coordinate system based on body in 0..200 (center 100), y in 0..250.
// We pack labelY tightly between 0 and 250 so the viewBox doesn't need to be artificially tall.
// This ensures the body takes up the full height of the SVG.
const organMap: Record<string, OrganCoordinate> = {
  // ── Right Side ──
  brain: { x: 50, y: 15, labelY: 10, label: "Brain", side: "right" },
  lung: { x: 60, y: 65, labelY: 38, label: "Lung", side: "right" },
  heart: { x: 50, y: 65, labelY: 66, label: "Heart", side: "right" },
  skin: { x: 95, y: 80, labelY: 94, label: "Skin", side: "right" },
  stomach: { x: 58, y: 85, labelY: 122, label: "Stomach", side: "right" },
  pancreas: { x: 55, y: 90, labelY: 150, label: "Pancreas", side: "right" },
  colon: { x: 50, y: 105, labelY: 178, label: "Colon / Intestine", side: "right" },
  endothelium: { x: 75, y: 130, labelY: 206, label: "Endothelium", side: "right" },
  muscle: { x: 75, y: 170, labelY: 234, label: "Muscle", side: "right" },

  // ── Left Side ──
  breast: { x: 38, y: 55, labelY: 10, label: "Breast", side: "left" }, 
  other_glandular: { x: 50, y: 35, labelY: 38, label: "Other Glandular", side: "left" },
  liver: { x: 42, y: 85, labelY: 66, label: "Liver", side: "left" },
  immune_blood: { x: 50, y: 48, labelY: 94, label: "Immune / Blood", side: "left" },
  kidney: { x: 40, y: 95, labelY: 122, label: "Kidney", side: "left" },
  ovary: { x: 40, y: 115, labelY: 150, label: "Female Reproductive", side: "left" },
  prostate: { x: 50, y: 125, labelY: 178, label: "Prostate", side: "left" },
  bone_marrow: { x: 25, y: 170, labelY: 206, label: "Bone Marrow", side: "left" },
  testis: { x: 50, y: 140, labelY: 234, label: "Male Reproductive", side: "left" },
};

export function normalizeTissueName(rawName: string): string | null {
  const name = rawName.toLowerCase();
  
  // Brain / CNS / Eye
  if (["brain", "cortex", "cerebellum", "ependymal", "astrocyte", "glial", "neuron", "microglia", "oligodendrocyte", "purkinje", "schwann", "meningeal", "pituitary", "hypothalamus", "somatotroph", "corticotroph", "gonadotroph", "thyrotroph", "pituicytes", "choroid plexus", "bergmann glia", "retinal", "cone", "rod", "ocular", "müller glia", "amacrine"].some(x => name.includes(x))) return "brain";
  // Breast
  if (["breast", "brca", "myoepithelial", "lactotroph"].some(x => name.includes(x))) return "breast";
  // Lung / Respiratory
  if (["lung", "alveolar", "pneumocyte", "bronchial", "respiratory", "club cells"].some(x => name.includes(x))) return "lung";
  // Heart
  if (["heart", "cardiac", "cardiomyocyte", "epicardial"].some(x => name.includes(x))) return "heart";
  // Liver
  if (["liver", "hepatic", "hepatocyte", "kupffer", "cholangiocyte", "ito cell", "gallbladder"].some(x => name.includes(x))) return "liver";
  // Stomach
  if (["stomach", "gastric", "parietal cell", "chief cell", "foveolar", "mucous neck"].some(x => name.includes(x))) return "stomach";
  // Kidney
  if (["kidney", "ccrcc", "renal", "podocyte", "mesangial", "macula densa", "proximal tubule", "distal tubule", "henle", "papillary tip", "distal convoluted"].some(x => name.includes(x))) return "kidney";
  // Colon / Intestine / Esophagus
  if (["colon", "colorectal", "intestine", "intestinal", "enterocyte", "goblet", "paneth", "duodenum", "jejunum", "ileum", "rectum", "tuft", "enteric", "submucosal glandular", "esophageal"].some(x => name.includes(x))) return "colon";
  // Pancreas
  if (["pancreas", "pancreatic", "islet", "beta cell", "alpha cell", "acinar"].some(x => name.includes(x))) return "pancreas";
  // Prostate
  if (["prostate", "prostatic"].some(x => name.includes(x))) return "prostate";
  // Ovary & Female Genital Tract
  if (["ovary", "ovarian", "follicle", "luteal", "fallopian", "endometrial", "cervix", "cervical", "vagina", "decidual", "oocyte", "granulosa"].some(x => name.includes(x))) return "ovary";
  
  // ── NEW SYSTEMS ──
  // Blood / Immune / Thymus
  if (["blood", "neutrophil", "macrophage", "monocyte", "t-cell", "b-cell", "nk-cell", "lymphocyte", "leukocyte", "erythrocyte", "platelet", "eosinophil", "basophil", "dendritic", "plasma cell", "thymocyte", "lymph node", "tonsil", "spleen", "splenic", "mast cell", "innate lymphoid", "pdcs", "cdc", "thymic"].some(x => name.includes(x))) return "immune_blood";
  // Bone Marrow
  if (["bone marrow", "myeloid", "erythroid", "megakaryocyte", "osteoclast", "osteoblast", "hematopoietic"].some(x => name.includes(x))) return "bone_marrow";
  // Skin & Connective
  if (["skin", "melanocyte", "keratinocyte", "fibroblast", "adipocyte", "epidermis", "dermis", "hair follicle", "sebaceous", "sweat gland"].some(x => name.includes(x))) return "skin";
  // Muscle / Skeleton
  if (["muscle", "myocyte", "skeletal", "smooth muscle", "myosatellite", "myonuclei", "fibro-adipogenic"].some(x => name.includes(x))) return "muscle";
  // Male Reproductive
  if (["testis", "testicular", "spermatid", "spermatocyte", "spermatogonia", "leydig", "sertoli", "epididymal", "seminal", "peritubular myoid"].some(x => name.includes(x))) return "testis";
  // Endothelium / Vessels
  if (["endothelial", "vein", "artery", "capillary", "pericyte"].some(x => name.includes(x))) return "endothelium";
  // Other Epithelial/Glandular/Placenta
  if (["thyroid", "parathyroid", "adrenal", "salivary", "urothelial", "bladder", "mesothelial", "trophoblast", "placenta", "amnion", "hofbauer"].some(x => name.includes(x))) return "other_glandular";

  return null;
}

interface TissueVisualizerProps {
  cohortTissues: string[]; 
  hpaTissues: string[]; // Unused for now, as user wants all HPA tissues as grey background map
  className?: string;
}

export default function TissueVisualizer({ cohortTissues, hpaTissues, className }: TissueVisualizerProps) {
  const cohortKeys = useMemo(() => {
    const keys = new Set<string>();
    cohortTissues.forEach((t) => {
      const norm = normalizeTissueName(t);
      if (norm) keys.add(norm);
    });
    return Array.from(keys);
  }, [cohortTissues]);

  const hpaKeys = useMemo(() => {
    const keys = new Set<string>();
    hpaTissues.forEach((t) => {
      const norm = normalizeTissueName(t);
      if (norm) keys.add(norm);
    });
    return Array.from(keys);
  }, [hpaTissues]);

  const allKeys = useMemo(() => {
    return Array.from(new Set([...cohortKeys, ...hpaKeys]));
  }, [cohortKeys, hpaKeys]);

  return (
    <div className={cn("relative flex items-center justify-center", className)}>
      <svg
        viewBox="-260 0 800 250"
        className="w-full h-auto max-h-[86px] opacity-90 overflow-visible"
      >
        {/* Minimalist Human Silhouette Path, made skinnier */}
        <path
          d="M100 5 C92 5, 87 12, 87 20 C87 28, 93 32, 93 32 L83 38 C75 42, 60 55, 55 70 C50 85, 45 110, 42 120 C40 128, 48 132, 52 125 C58 115, 65 95, 68 85 L73 130 C75 140, 72 170, 70 200 C68 220, 65 240, 75 240 C85 240, 85 220, 88 190 C90 160, 95 140, 100 140 C105 140, 110 160, 112 190 C115 220, 115 240, 125 240 C135 240, 132 220, 130 200 C128 170, 125 140, 127 130 L132 85 C135 95, 142 115, 148 125 C152 132, 160 128, 158 120 C155 110, 150 85, 145 70 C140 55, 125 42, 117 38 L107 32 C107 32, 113 28, 113 20 C113 12, 108 5, 100 5 Z"
          fill="currentColor"
          className="text-muted-foreground/15"
          stroke="currentColor"
          strokeWidth="3.5"
          strokeLinejoin="round"
          style={{ transformOrigin: "100px 125px", transform: "scaleX(0.70)" }}
        />

        {/* Render active organs */}
        {allKeys.map((key) => {
          const org = organMap[key];
          if (!org) return null;
          
          const isCohort = cohortKeys.includes(key);

          // x=100 is the center, offset mapping slightly
          const cx = org.x + 50; 
          const cy = org.y;
          const isLeft = org.side === "left";
          
          // Labels pushed far horizontally, safely inside the -260 / 540 boundaries
          const labelX = isLeft ? cx - 150 : cx + 110;
          const labelY = org.labelY;

          // Polyline points for a sleek elbow connector
          const elbowX = isLeft ? cx - 30 : cx + 30;
          const elbowX2 = isLeft ? labelX + 15 : labelX - 15;
          const points = `${cx},${cy} ${elbowX},${cy} ${elbowX2},${labelY} ${labelX},${labelY}`;

          // Styles based on cohort vs HPA
          const dotFillClass = isCohort ? "fill-primary" : "fill-muted-foreground/30";
          const lineStrokeClass = isCohort ? "text-primary/70" : "text-muted-foreground/25";
          const textClass = isCohort ? "fill-primary font-bold" : "fill-muted-foreground/40 font-medium";

          return (
            <g key={key} className="animate-in fade-in zoom-in duration-500">
              {/* Callout Line */}
              <polyline
                points={points}
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                className={lineStrokeClass}
                strokeDasharray={isCohort ? "none" : "5 5"}
              />
              
              {/* Pulsing Dot (Cohort only) */}
              {isCohort && (
                <circle
                  cx={cx}
                  cy={cy}
                  r="14"
                  className="fill-primary/30 animate-ping origin-center"
                  style={{ transformOrigin: `${cx}px ${cy}px`, animationDuration: '2s' }}
                />
              )}
              
              {/* Solid Dot */}
              <circle
                cx={cx}
                cy={cy}
                r="7"
                className={`${dotFillClass} stroke-background`}
                strokeWidth="2.5"
              />
              
              {/* For breast, we also want a dot on the other side */}
              {key === 'breast' && (
                <>
                  {isCohort && (
                    <circle
                      cx={cx + 24}
                      cy={cy}
                      r="14"
                      className="fill-primary/30 animate-ping origin-center"
                      style={{ transformOrigin: `${cx + 24}px ${cy}px`, animationDuration: '2s' }}
                    />
                  )}
                  <circle
                    cx={cx + 24}
                    cy={cy}
                    r="7"
                    className={`${dotFillClass} stroke-background`}
                    strokeWidth="2.5"
                  />
                </>
              )}

              {/* Text Label */}
              <text
                x={labelX + (isLeft ? -10 : 10)}
                y={labelY + 11} // Adjusted vertical alignment for font size
                textAnchor={isLeft ? "end" : "start"}
                fontSize="32"
                className={`${textClass} tracking-widest uppercase`}
              >
                {org.label}
              </text>
            </g>
          );
        })}
      </svg>
      {/* Legend */}
      <div className="absolute bottom-[-16px] left-0 right-0 flex justify-center gap-6 text-[10px] text-muted-foreground/80 font-medium">
        <div className="flex items-center gap-1.5">
          <span className="w-2.5 h-2.5 rounded-full bg-primary inline-block" />
          <span>Datasets</span>
        </div>
        <div className="flex items-center gap-1.5">
          <span className="w-2.5 h-2.5 rounded-full bg-muted-foreground/30 inline-block" />
          <span>HPA Reference</span>
        </div>
      </div>
    </div>
  );
}
