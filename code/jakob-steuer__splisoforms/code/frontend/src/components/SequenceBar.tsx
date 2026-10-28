"use client";

import React, { useMemo, useState, useEffect } from "react";
import { Copy, Download, Check, Info } from "lucide-react";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";

export interface JunctionEvent {
  id: number;
  junction_type: string;
  transcript_pos?: number;
  is_novel: boolean;
}

interface SequenceBarProps {
  sequence: string;
  junctionEvents?: JunctionEvent[];
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onClick?: (resNo: number) => void;
  onSelectionChange?: (residues: number[]) => void;
}

const JUNCTION_STYLES: Record<string, string> = {
  novel_exon: "bg-blue-500/30 ring-1 ring-blue-500/50",
  skipped_exon: "opacity-40 line-through decoration-slate-500",
  alt_5_splice: "bg-yellow-500/40 ring-1 ring-yellow-500/60",
  alt_3_splice: "bg-yellow-500/40 ring-1 ring-yellow-500/60",
  intron_retention: "border-b-2 border-dotted border-purple-500",
};

export function SequenceBar({ 
  sequence, 
  junctionEvents = [], 
  hoveredResidues = [], 
  selectedResidues = [], 
  onHover, 
  onClick, 
  onSelectionChange 
}: SequenceBarProps) {
  const [selectionStart, setSelectionStart] = useState<number | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [copied, setCopied] = useState(false);

  // Global mouse up to end dragging
  useEffect(() => {
    const handleGlobalMouseUp = () => setIsDragging(false);
    window.addEventListener("mouseup", handleGlobalMouseUp);
    return () => window.removeEventListener("mouseup", handleGlobalMouseUp);
  }, []);

  const handleCopy = () => {
    navigator.clipboard.writeText(sequence);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handleDownload = () => {
    const blob = new Blob([`>Sequence\n${sequence}`], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "sequence.fasta";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const handleResidueMouseDown = (e: React.MouseEvent, resNo: number) => {
    setIsDragging(true);
    setSelectionStart(resNo);

    if (e.shiftKey && selectionStart !== null) {
      const start = Math.min(selectionStart, resNo);
      const end = Math.max(selectionStart, resNo);
      const range = [];
      for (let i = start; i <= end; i++) range.push(i);
      
      if (e.metaKey || e.ctrlKey) {
         const combined = Array.from(new Set([...selectedResidues, ...range]));
         if (onSelectionChange) onSelectionChange(combined);
      } else {
         if (onSelectionChange) onSelectionChange(range);
      }
    } else if (e.metaKey || e.ctrlKey) {
      if (onSelectionChange) {
         if (selectedResidues.includes(resNo)) {
             onSelectionChange(selectedResidues.filter(r => r !== resNo));
         } else {
             onSelectionChange([...selectedResidues, resNo]);
         }
      }
    } else {
      if (onClick) onClick(resNo);
      if (onSelectionChange) onSelectionChange([resNo]);
    }
  };

  const handleResidueMouseEnter = (resNo: number) => {
    if (onHover) onHover([resNo]);
    
    if (isDragging && selectionStart !== null) {
      const start = Math.min(selectionStart, resNo);
      const end = Math.max(selectionStart, resNo);
      const range = [];
      for (let i = start; i <= end; i++) range.push(i);
      if (onSelectionChange) onSelectionChange(range);
    }
  };

  const CHUNK_SIZE = 10;
  const chunks = [];
  for (let i = 0; i < sequence.length; i += CHUNK_SIZE) {
    chunks.push({
      start: i + 1,
      sequence: sequence.slice(i, i + CHUNK_SIZE)
    });
  }

  return (
    <div className="w-full bg-card border border-border rounded-lg p-4 shadow-inner">
      <div className="flex justify-end items-center mb-4 pb-2 border-b border-border/50">
        <div className="flex gap-2">
          <button 
            onClick={handleCopy} 
            className="flex items-center gap-1.5 px-2 py-1 bg-muted hover:bg-muted/80 text-muted-foreground text-xs rounded transition-colors"
          >
            {copied ? <Check className="w-3 h-3 text-emerald-400" /> : <Copy className="w-3 h-3" />}
            {copied ? "Copied" : "Copy"}
          </button>
          <button 
            onClick={handleDownload} 
            className="flex items-center gap-1.5 px-2 py-1 bg-primary/20 hover:bg-primary/30 text-primary text-xs rounded transition-colors"
          >
            <Download className="w-3 h-3" />
            Download FASTA
          </button>
        </div>
      </div>

      <TooltipProvider>
        <div 
          className="flex flex-wrap gap-x-4 gap-y-6 font-mono text-sm leading-none max-h-[400px] overflow-y-auto pr-2 custom-scrollbar select-none"
          onMouseLeave={() => {
            setIsDragging(false);
            if (onHover) onHover([]);
          }}
        >
          {chunks.map((chunk, cIdx) => (
            <div key={cIdx} className="flex flex-col gap-1.5">
              {/* Numbering directly above each block */}
              <div className="text-[10px] text-muted-foreground/60 shrink-0 select-none tracking-widest pl-0.5">
                {chunk.start}
              </div>
              
              {/* Sequence Block of 10 */}
              <div className="flex">
                {chunk.sequence.split('').map((char, charIdx) => {
                  const resNo = chunk.start + charIdx;
                  const eventsAtPos = junctionEvents.filter(e => e.transcript_pos === resNo);
                  const isHovered = hoveredResidues?.includes(resNo) ?? false;
                  const isSelected = selectedResidues?.includes(resNo) ?? false;
                  
                  const junctionClass = eventsAtPos.length > 0 ? JUNCTION_STYLES[eventsAtPos[0].junction_type] || "bg-primary/10" : "";

                  const item = (
                    <span
                      key={charIdx}
                      className={`cursor-crosshair inline-flex items-center justify-center w-[1.2ch] h-6 transition-all duration-150 ${
                        isSelected
                          ? "bg-primary text-primary-foreground font-bold rounded scale-110 shadow-lg z-10"
                          : isHovered
                          ? "bg-secondary/40 text-secondary-foreground font-bold rounded scale-110"
                          : junctionClass
                          ? `${junctionClass} rounded`
                          : "text-muted-foreground hover:text-foreground hover:bg-muted"
                      }`}
                      onMouseEnter={() => handleResidueMouseEnter(resNo)}
                      onMouseDown={(e) => handleResidueMouseDown(e, resNo)}
                    >
                      {char}
                    </span>
                  );

                  if (eventsAtPos.length > 0) {
                    return (
                      <Tooltip key={charIdx}>
                        <TooltipTrigger>
                          {item}
                        </TooltipTrigger>
                        <TooltipContent className="bg-popover border-border shadow-xl p-2 text-xs">
                          <div className="space-y-1">
                            <div className="font-bold flex items-center gap-1 text-primary">
                              <Info className="w-3 h-3" /> Junction Detail
                            </div>
                            {eventsAtPos.map((e, idx) => (
                              <div key={idx} className="text-foreground">
                                <span className="capitalize">{e.junction_type.replace('_', ' ')}</span>
                                {e.is_novel && <span className="ml-1 text-[10px] bg-red-500/20 text-red-400 px-1 rounded">Novel</span>}
                              </div>
                            ))}
                            <div className="text-muted-foreground text-[10px]">Position: {resNo}</div>
                          </div>
                        </TooltipContent>
                      </Tooltip>
                    );
                  }

                  return item;
                })}
              </div>
            </div>
          ))}
        </div>
      </TooltipProvider>
    </div>
  );
}

