"use client";

import { useState, useEffect } from "react";
import { getApiUrl } from "@/lib/api";

export function StatsBanner() {
  const [stats, setStats] = useState({ isoform_count: 0, gene_count: 0 });

  useEffect(() => {
    const fetchStats = async () => {
      try {
        const apiUrl = getApiUrl();
        const res = await fetch(`${apiUrl}/api/stats`);
        if (res.ok) {
          const data = await res.json();
          setStats(data);
        }
      } catch (err) {
        console.error("Failed to fetch stats", err);
      }
    };
    fetchStats();
    
    const interval = setInterval(fetchStats, 10000);
    return () => clearInterval(interval);
  }, []);

  if (stats.isoform_count === 0) return null;

  return (
    <div className="mb-12 text-md font-semibold text-primary animate-pulse flex items-center justify-center gap-3 bg-background/80 border border-primary/40 px-6 py-3 rounded-full shadow-[0_0_20px_rgba(var(--primary-rgb),0.15)] backdrop-blur-sm">
      <span className="relative flex h-3 w-3">
        <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-primary opacity-75"></span>
        <span className="relative inline-flex rounded-full h-3 w-3 bg-primary"></span>
      </span>
      Analyzing {stats.isoform_count.toLocaleString()} structurally modelled isoforms across {stats.gene_count.toLocaleString()} genes — reference, tumour tissue, and organoids.
    </div>
  );
}
