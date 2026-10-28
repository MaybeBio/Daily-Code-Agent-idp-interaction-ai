"use client";

import { useEffect, useId, useRef, useState } from "react";

/**
 * 2D skeletal formula rendered from SMILES, client-side.
 *
 * Drawn locally with smiles-drawer — deliberately NOT PubChem's depiction image, which
 * would be an external request per drug on every page view of a public site.
 *
 * The library mutates the <svg> in place, so it needs a ref and a real DOM node; it also
 * fails softly on unparseable SMILES (calls onError rather than throwing), which we surface
 * as an explicit message instead of an empty box.
 */
export default function MoleculeSkeleton({
  smiles,
  className = "",
}: {
  smiles?: string | null;
  className?: string;
}) {
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [error, setError] = useState<string | null>(null);
  // Unique id per instance: smiles-drawer keys its <defs> off the target element.
  const uid = useId().replace(/:/g, "");

  useEffect(() => {
    const svg = svgRef.current;
    if (!svg || !smiles) return;

    let cancelled = false;
    setError(null);

    // Dynamic import: the library touches the DOM at module scope, so it must not run
    // during SSR.
    import("smiles-drawer")
      .then(({ default: SmilesDrawer }) => {
        if (cancelled || !svgRef.current) return;

        const drawer = new SmilesDrawer.SvgDrawer({
          width: 340,
          height: 240,
          padding: 12,
          bondThickness: 1.1,
          terminalCarbons: true,
          themes: {
            // Muted, theme-agnostic palette: heteroatoms keep their conventional colours so
            // the diagram stays readable, but carbon/bonds follow the page's foreground.
            splisoforms: {
              C: "#64748b",
              N: "#3b82f6",
              O: "#ef4444",
              S: "#eab308",
              P: "#f97316",
              F: "#22c55e",
              CL: "#22c55e",
              BR: "#a16207",
              I: "#a855f7",
              BACKGROUND: "#00000000",
            },
          },
        });

        SmilesDrawer.parse(
          smiles,
          (tree) => {
            if (cancelled || !svgRef.current) return;
            try {
              drawer.draw(tree, svgRef.current, "splisoforms");
            } catch {
              setError("Could not render this structure.");
            }
          },
          () => {
            if (!cancelled) setError("Unparseable SMILES.");
          },
        );
      })
      .catch(() => {
        if (!cancelled) setError("Structure renderer failed to load.");
      });

    return () => {
      cancelled = true;
    };
  }, [smiles, uid]);

  if (!smiles) {
    return (
      <div className={`flex items-center justify-center text-xs text-muted-foreground ${className}`}>
        No structure available
      </div>
    );
  }

  return (
    <div className={`relative flex items-center justify-center ${className}`}>
      <svg
        ref={svgRef}
        id={`sk-${uid}`}
        className="w-full h-full"
        aria-label="2D skeletal formula"
      />
      {error && (
        <span className="absolute text-xs text-muted-foreground">{error}</span>
      )}
    </div>
  );
}
