"use client";

import dynamic from "next/dynamic";
import { useEffect, useState } from "react";
import { useTheme } from "@/context/ThemeContext";
import { getApiUrl } from "@/lib/api";

// Dynamically load Plotly to prevent Next.js SSR issues with browser APIs.
const Plot = dynamic(() => import("react-plotly.js"), { 
  ssr: false, 
  loading: () => (
    <div className="w-full h-full min-h-[400px] flex items-center justify-center bg-muted/50 rounded-lg border border-border animate-pulse">
      <span className="text-muted-foreground font-medium">Loading PAE Matrix...</span>
    </div>
  ) 
});

export default function PAEMatrix({ 
  paeUrl, 
  sequenceLength,
  onHover 
}: { 
  paeUrl?: string; 
  sequenceLength?: number;
  onHover?: (resPairs: number[]) => void;
}) {
  const [data, setData] = useState<{z: number[][]} | null>(null);
  const { theme } = useTheme();
  const isDark = theme === "dark";
  const gridColor = isDark ? "#334155" : "#e2e8f0";
  const textColor = isDark ? "#94A3B8" : "#475569";

  useEffect(() => {
    const fetchData = async () => {
      if (paeUrl) {
        const BASE_URL = getApiUrl().replace(/\/api$/, '');
        
        const joinUrl = (base: string, path: string) => {
          if (path.startsWith('http')) return path;
          const cleanBase = base.replace(/\/+$/, '');
          const cleanPath = path.replace(/^\/+/, '');
          return `${cleanBase}/${cleanPath}`;
        };
        
        const finalUrl = joinUrl(BASE_URL, paeUrl);
        
        try {
           const res = await fetch(finalUrl);
           if (res.ok) {
              const json = await res.json();
              let paeData = null;
              
              if (Array.isArray(json) && json[0] && json[0].predicted_aligned_error) {
                 paeData = json[0].predicted_aligned_error;
              } else if (json.predicted_aligned_error) {
                 paeData = json.predicted_aligned_error;
              } else if (json.pae) {
                 paeData = json.pae;
              } else if (Array.isArray(json) && Array.isArray(json[0])) {
                 paeData = json;
              }
              
              if (paeData) {
                 setData({ z: paeData });
                 return;
              }
           }
        } catch(e) {
           console.error("Failed to load PAE:", e);
        }
      }
      
      // Fallback mock ONLY if fetch fails or no URL
      const size = 50;
      const mockZ = Array.from({length: size}, (_, i) => 
        Array.from({length: size}, (_, j) => {
           if (i < 20 && j < 20) return Math.random() * 5 + 2; 
           if (i > 30 && j > 30) return Math.random() * 5 + 2;
           return Math.random() * 20 + 10;
        })
      );
      setData({ z: mockZ });
    };
    fetchData();
  }, [paeUrl]);

  const [sizeWarning, setSizeWarning] = useState(false);
  useEffect(() => {
    if (data && sequenceLength && data.z.length !== sequenceLength && data.z.length > 50) {
      // 50 is our mock length. If real data length conflicts with seq length, warn
      setSizeWarning(true);
    }
  }, [data, sequenceLength]);

  if (!data) return null;

  return (
    <div className="w-full h-full min-h-[400px] bg-card/30 rounded-lg overflow-hidden border border-border relative shadow-inner">
      {sizeWarning && (
        <div className="absolute top-2 left-2 z-10 bg-destructive/80 text-white text-xs px-2 py-1 rounded shadow-lg backdrop-blur">
          Warning: PAE size mismatch with sequence length.
        </div>
      )}
      <Plot
        onHover={(e) => {
           if (onHover && e.points && e.points[0]) {
              onHover([e.points[0].x as number, e.points[0].y as number]);
           }
        }}
        onUnhover={() => {
           if (onHover) onHover([]);
        }}
        data={[
          {
            z: data.z,
            type: "heatmap",
            colorscale: [
              [0, "#004D40"], // Deep Green (low error)
              [1, "#FFFFFF"]  // White (high error)
            ],
            showscale: true,
            colorbar: {
               title: { 
                 text: "Expected Error (Å)",
                 font: { color: textColor }
               },
               tickfont: { color: textColor },
            }
          }
        ]}
        layout={{
          autosize: true,
          margin: { t: 20, l: 70, r: 90, b: 70 },
          paper_bgcolor: "transparent",
          plot_bgcolor: "transparent",
          font: { color: textColor },
          xaxis: {
             title: { text: "Scored Residue", standoff: 12 },
             gridcolor: gridColor,
             showspikes: true,
             spikemode: "across",
             spikedash: "solid",
             spikecolor: "#10B981",
             spikethickness: 1
          },
          yaxis: {
             title: { text: "Aligned Residue", standoff: 12 },
             gridcolor: gridColor,
             autorange: "reversed",
             // PAE is an N×N residue matrix — lock the data aspect ratio to 1:1 so
             // the heatmap always renders square regardless of container width.
             scaleanchor: "x",
             scaleratio: 1,
             constrain: "domain",
             showspikes: true,
             spikemode: "across",
             spikedash: "solid",
             spikecolor: "#10B981",
             spikethickness: 1
          },
          hovermode: "closest",
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%", position: "absolute", inset: 0 }}
      />
    </div>
  );
}
