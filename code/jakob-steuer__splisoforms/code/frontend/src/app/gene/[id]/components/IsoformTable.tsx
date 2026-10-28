"use client";

import { useState, useCallback } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { ArrowRight, Box, ShieldAlert, Flag, Filter, ChevronUp, ChevronDown, ChevronsUpDown, CheckSquare, Square } from "lucide-react";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { type SortKey, type SortOrder } from "./GeneDashboardClient";
import { getApiUrl } from "@/lib/api";

// Comparison view caps at 5 same-gene isoforms (overlay readability).
// Selection itself is uncapped — exceeding the cap just informs the Compare
// button to show a "deselect N to compare" hint instead of becoming active.
const MAX_COMPARE = 5;

export default function IsoformTable({
  isoforms: initialIsoforms,
  selectedIds = [],
  onToggleSelect,
  onSetSelection,
  sortKey,
  sortOrder,
  onSortChange,
}: {
  isoforms: any[],
  selectedIds?: number[],
  onToggleSelect?: (id: number) => void,
  /** Replace the entire selection at once — used for Select All / Deselect All */
  onSetSelection?: (ids: number[]) => void,
  sortKey: SortKey,
  sortOrder: SortOrder,
  onSortChange: (key: SortKey, order: SortOrder) => void,
}) {
  const [filter, setFilter] = useState("");
  const [sourceFilter, setSourceFilter] = useState<string>("all");
  const [statusFilter, setStatusFilter] = useState<string>("all");
  const [isoforms, setIsoforms] = useState(initialIsoforms);

  // All unique source names across all isoforms, including multi-source memberships
  const sources = Array.from(new Set(
    initialIsoforms.flatMap((i) => {
      const memberships: any[] = i.dataset_presence ?? [];
      if (memberships.length > 0) return memberships.map((m: any) => m.source_name as string);
      return [i.source?.name || "Unknown"];
    })
  )) as string[];

  const handleSort = (key: SortKey) => {
    if (sortKey === key) {
      onSortChange(key, sortOrder === "asc" ? "desc" : "asc");
    } else {
      onSortChange(key, "asc");
    }
  };

  const SortIcon = ({ k }: { k: SortKey }) => {
    if (sortKey !== k) return <ChevronsUpDown className="w-3 h-3 ml-1.5 opacity-30" />;
    return sortOrder === "asc" 
      ? <ChevronUp className="w-3 h-3 ml-1.5 text-primary" /> 
      : <ChevronDown className="w-3 h-3 ml-1.5 text-primary" />;
  };

  const toggleFlag = useCallback(async (id: number) => {
    // ... same as before
    setIsoforms(prev =>
      prev.map(iso => iso.id === id ? { ...iso, is_flagged: !iso.is_flagged } : iso)
    );
    try {
      const apiUrl = getApiUrl();
      await fetch(`${apiUrl}/api/isoform/${id}/flag`, { method: "PATCH" });
    } catch {
      setIsoforms(prev =>
        prev.map(iso => iso.id === id ? { ...iso, is_flagged: !iso.is_flagged } : iso)
      );
    }
  }, []);

  if (!isoforms || isoforms.length === 0) {
    return (
      <div className="bg-card/40 border border-primary/20 rounded-2xl p-10 flex flex-col items-center justify-center text-center">
        <Box className="w-12 h-12 text-muted-foreground mb-4 opacity-50" />
        <h3 className="text-xl font-medium text-foreground mb-2">No Transcripts Found</h3>
        <p className="text-muted-foreground max-w-sm">
          There are currently no mapped isoforms for this gene in the database.
        </p>
      </div>
    );
  }

  const filteredIsoforms = isoforms
    .filter((iso) => {
      const q = filter.toLowerCase();
      const matchesText = iso.transcript_id.toLowerCase().includes(q) ||
        (iso.transcript_name && iso.transcript_name.toLowerCase().includes(q)) ||
        (iso.metadata_record?.description && iso.metadata_record.description.toLowerCase().includes(q));
      const isoSources: string[] = (iso.dataset_presence ?? []).length > 0
        ? (iso.dataset_presence as any[]).map((m: any) => m.source_name as string)
        : [iso.source?.name || "Unknown"];
      const matchesSource = sourceFilter === "all" || isoSources.includes(sourceFilter);
      const matchesStatus = statusFilter === "all" || 
        (statusFilter === "canonical" && iso.is_canonical) || 
        (statusFilter === "novel" && !iso.is_canonical);
      return matchesText && matchesSource && matchesStatus;
    })
    .sort((a, b) => {
       // 1. Canonical always pinned to top
       if (a.is_canonical !== b.is_canonical) return a.is_canonical ? -1 : 1;
       
       // 2. Perform dynamic sort
       let comparison = 0;
       if (sortKey === "id") {
         comparison = a.transcript_id.localeCompare(b.transcript_id);
       } else if (sortKey === "name") {
         // Unnamed novel rows (no close GENCODE match) sort last rather than first.
         const nameA = a.transcript_name ?? "";
         const nameB = b.transcript_name ?? "";
         if (!nameA !== !nameB) return nameA ? -1 : 1;
         comparison = nameA.localeCompare(nameB, undefined, { numeric: true });
       } else if (sortKey === "length") {
         const lenA = a.sequence?.length || 0;
         const lenB = b.sequence?.length || 0;
         comparison = lenA - lenB;
       } else if (sortKey === "source") {
         const srcA = a.source?.name || "Unknown";
         const srcB = b.source?.name || "Unknown";
         comparison = srcA.localeCompare(srcB);
       }

       if (comparison !== 0) {
         return sortOrder === "asc" ? comparison : -comparison;
       }

       // 3. Fallback to flag then ID for stability
       if (a.is_flagged !== b.is_flagged) return a.is_flagged ? -1 : 1;
       return a.transcript_id.localeCompare(b.transcript_id);
    });

  // Derived selection helpers
  const filteredIds = filteredIsoforms.map(iso => iso.id);
  const allFilteredSelected = filteredIds.length > 0 && filteredIds.every(id => selectedIds.includes(id));

  const handleSelectToggleAll = () => {
    if (!onSetSelection) return;
    if (allFilteredSelected) {
      // Deselect all currently filtered; keep any selections outside the filter
      onSetSelection(selectedIds.filter(id => !filteredIds.includes(id)));
    } else {
      // Add all filtered ids that aren't already selected (no cap — Compare button
      // alone surfaces the 5-isoform overlay limit).
      onSetSelection(Array.from(new Set([...selectedIds, ...filteredIds])));
    }
  };

  // Derived Compare-button state.
  const selectionCount = selectedIds.length;
  const compareOverLimit = selectionCount > MAX_COMPARE;
  const compareEnabled = selectionCount >= 2 && !compareOverLimit;

  return (
    <div className="bg-card/60 backdrop-blur-md border border-primary/20 shadow-xl rounded-2xl flex flex-col overflow-hidden">
      <div className="p-6 border-b border-border/50 flex flex-col gap-4">
        {/* Row 1: title block (left) + filter cluster (right) */}
        <div className="flex flex-col lg:flex-row justify-between items-start lg:items-center gap-4">
          <div className="space-y-1">
            <h2 className="text-xl font-semibold text-foreground">Mapped Transcripts</h2>
            <p className="text-sm text-muted-foreground">
              {isoforms.length} isoform{isoforms.length === 1 ? "" : "s"} for this gene · pick 2–{MAX_COMPARE} to compare side-by-side.
            </p>
          </div>
          <div className="flex gap-2 flex-wrap items-center">
            {sources.length > 1 && (
              <select
                className="bg-background/50 border border-border rounded-md px-3 py-1.5 text-sm text-muted-foreground focus:outline-none"
                value={sourceFilter}
                onChange={(e) => setSourceFilter(e.target.value)}
              >
                <option value="all">All Sources</option>
                {sources.map(s => <option key={s} value={s}>{s}</option>)}
              </select>
            )}
            <select
              className="bg-background/50 border border-border rounded-md px-3 py-1.5 text-sm text-muted-foreground focus:outline-none"
              value={statusFilter}
              onChange={(e) => setStatusFilter(e.target.value)}
            >
              <option value="all">All Status</option>
              <option value="canonical">Canonical Only</option>
              <option value="novel">Novel Only</option>
            </select>
            <Input
              placeholder="Filter by name or ID..."
              className="bg-background/50 border-border/50 w-48"
              value={filter}
              onChange={(e) => setFilter(e.target.value)}
            />
          </div>
        </div>

        {/* Row 2: selection toolbar — visible only when something is selected,
            keeps the Compare action visually grouped with the select-all toggle
            and the selection count, separated from the filter cluster. */}
        {(onSetSelection || selectionCount > 0) && (
          <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3 pt-3 border-t border-border/30">
            <div className="flex items-center gap-3 flex-wrap">
              <span className="text-xs uppercase tracking-wider text-muted-foreground">
                <span className="font-semibold text-foreground">{selectionCount}</span>
                {" "}selected
              </span>
              {onSetSelection && filteredIds.length > 0 && (
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={handleSelectToggleAll}
                  className="h-8 gap-1.5 text-xs text-muted-foreground hover:text-primary hover:bg-primary/10 border border-border/60"
                >
                  {allFilteredSelected
                    ? <><Square className="w-3.5 h-3.5" /> Deselect All</>
                    : <><CheckSquare className="w-3.5 h-3.5" /> Select All</>
                  }
                </Button>
              )}
              {onSetSelection && selectionCount > 0 && (
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={() => onSetSelection([])}
                  className="h-8 px-2 text-xs text-muted-foreground hover:text-rose-400 hover:bg-rose-500/10"
                >
                  Clear
                </Button>
              )}
            </div>

            {/* Compare button — three visual states */}
            <div className="flex items-center gap-2">
              {compareEnabled ? (
                <Link href={`/compare?ids=${selectedIds.join(",")}`}>
                  <Button
                    variant="secondary"
                    size="sm"
                    className="bg-primary/20 hover:bg-primary/30 text-primary border-primary/30 h-8 gap-2"
                  >
                    <Filter className="w-3.5 h-3.5" />
                    Compare Selection ({selectionCount})
                  </Button>
                </Link>
              ) : compareOverLimit ? (
                <Button
                  variant="secondary"
                  size="sm"
                  disabled
                  className="h-8 gap-2 bg-amber-500/10 text-amber-500 border-amber-500/30 cursor-not-allowed"
                  title={`Comparison overlay is capped at ${MAX_COMPARE} isoforms. Deselect ${selectionCount - MAX_COMPARE} to enable.`}
                >
                  <Filter className="w-3.5 h-3.5" />
                  Compare limited to {MAX_COMPARE} · deselect {selectionCount - MAX_COMPARE}
                </Button>
              ) : (
                <span className="text-xs text-muted-foreground italic">
                  Select 2 or more to compare
                </span>
              )}
            </div>
          </div>
        )}
      </div>
      
      <div className="overflow-x-auto custom-scrollbar">
        <table className="w-full text-left text-sm text-muted-foreground whitespace-nowrap">
          <thead className="text-xs uppercase bg-muted/30 text-foreground/70 tracking-wider">
            <tr>
              <th scope="col" className="px-4 py-4 font-medium w-10">Flag</th>
              <th scope="col" className="px-4 py-4 font-medium w-10">Overlay</th>
              
              <th onClick={() => handleSort("name")} scope="col" className="px-6 py-4 font-medium cursor-pointer hover:text-foreground transition-colors">
                <div className="flex items-center">Transcript Name <SortIcon k="name" /></div>
              </th>

              <th onClick={() => handleSort("id")} scope="col" className="px-6 py-4 font-medium cursor-pointer hover:text-foreground transition-colors">
                <div className="flex items-center">Transcript ID <SortIcon k="id" /></div>
              </th>

              <th onClick={() => handleSort("length")} scope="col" className="px-6 py-4 font-medium cursor-pointer hover:text-foreground transition-colors">
                <div className="flex items-center">Length (aa) <SortIcon k="length" /></div>
              </th>
              
              <th onClick={() => handleSort("source")} scope="col" className="px-6 py-4 font-medium cursor-pointer hover:text-foreground transition-colors">
                <div className="flex items-center">Source <SortIcon k="source" /></div>
              </th>
              
              <th scope="col" className="px-6 py-4 font-medium text-right">Action</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-border/30">
            {filteredIsoforms.map((iso) => {
              const memberships: any[] = iso.dataset_presence ?? [];
              const sourceNames: string[] = memberships.length > 0
                ? memberships.map((m: any) => m.source_name as string)
                : [iso.source?.name || "Unknown"];
              const isNMD = iso.nmd_prediction?.is_nmd;
              const isSelected = selectedIds.includes(iso.id);
              const isFlagged = iso.is_flagged;
              
              return (
                <tr key={iso.id} className={`hover:bg-muted/50 transition-colors group ${isFlagged ? 'bg-amber-500/5' : ''}`}>
                  <td className="px-4 py-4">
                    <button
                      onClick={() => toggleFlag(iso.id)}
                      title={isFlagged ? "Unmark" : "Flag as isoform of interest"}
                      className={`p-1.5 rounded transition-colors ${
                        isFlagged 
                          ? 'text-amber-400 bg-amber-500/15 hover:bg-amber-500/25' 
                          : 'text-slate-600 hover:text-amber-400 hover:bg-amber-500/10'
                      }`}
                    >
                      <Flag className="w-3.5 h-3.5" fill={isFlagged ? 'currentColor' : 'none'} />
                    </button>
                  </td>
                  <td className="px-4 py-4">
                    <input
                      type="checkbox"
                      className="w-4 h-4 rounded border-border bg-background/50 text-primary focus:ring-primary/50 cursor-pointer"
                      checked={isSelected}
                      onChange={() => onToggleSelect?.(iso.id)}
                    />
                  </td>
                  {/* Name is the readable handle (MET-201 / MET-206-NOV-a6d00a) and carries
                      the link; transcript_id stays visible as the stable source identifier.
                      Novel isoforms with no close GENCODE match have no name — fall back. */}
                  <td className="px-6 py-4 font-mono font-medium text-primary/90">
                    <Link href={`/isoform/${iso.id}`} className="hover:text-primary hover:underline transition-colors">
                      {iso.transcript_name || iso.transcript_id}
                    </Link>
                    {iso.is_canonical && (
                      <Badge variant="outline" className="ml-3 text-[10px] uppercase bg-amber-500/10 text-amber-500 border-amber-500/20">
                        Canonical
                      </Badge>
                    )}
                  </td>
                  <td className="px-6 py-4 font-mono text-xs text-muted-foreground whitespace-nowrap">
                    {iso.transcript_id}
                  </td>
                  <td className="px-6 py-4">
                    <div className="flex items-center gap-2">
                      <Badge variant="secondary" className="bg-muted text-muted-foreground">
                        {iso.sequence?.length || 0} aa
                      </Badge>
                      {iso.delta_length !== null && iso.delta_length !== undefined && !iso.is_canonical && (
                        <span className={`text-xs ${iso.delta_length > 0 ? "text-emerald-400" : iso.delta_length < 0 ? "text-rose-400" : "text-muted-foreground font-medium"}`}>
                          {iso.delta_length > 0 ? "+" : iso.delta_length === 0 ? "±" : ""}{iso.delta_length} aa
                        </span>
                      )}
                      {isNMD && (
                         <div className="flex items-center text-rose-400 bg-rose-500/10 px-2 py-0.5 rounded text-xs gap-1">
                           <ShieldAlert className="w-3 h-3" /> NMD Target
                         </div>
                      )}
                    </div>
                  </td>
                  <td className="px-6 py-4">
                    <div className="flex flex-wrap gap-1">
                      {sourceNames.map((name) => {
                        const shortName = name.replace(" Nanopore ESPRESSO", "").replace(" PacBio MAS-seq", "");
                        return (
                          <Link 
                            key={name} 
                            href="/documentation"
                            title={name}
                            className="text-xs uppercase tracking-wider px-2 py-1 rounded bg-indigo-500/10 text-indigo-400 border border-indigo-500/20 whitespace-nowrap hover:bg-indigo-500/20 transition-colors"
                          >
                            {shortName}
                          </Link>
                        );
                      })}
                    </div>
                  </td>
                  <td className="px-6 py-4 text-right">
                    <Link 
                      href={`/isoform/${iso.id}`}
                      className="inline-flex items-center justify-center size-8 rounded-full bg-primary/10 text-primary hover:bg-primary hover:text-primary-foreground transition-all duration-300 transform group-hover:-translate-x-1"
                    >
                      <ArrowRight className="w-4 h-4" />
                    </Link>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
        {filteredIsoforms.length === 0 && filter !== "" && (
           <div className="p-8 text-center text-muted-foreground">
              No transcripts match your search term "{filter}".
           </div>
        )}
      </div>
    </div>
  );
}
