"use client";

import dynamic from "next/dynamic";
import { ExternalLink, FlaskConical } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import MoleculeSkeleton from "@/components/MoleculeSkeleton";

const StructureViewer = dynamic(() => import("@/components/StructureViewer"), { ssr: false });

export interface DrugSummary {
  id: number;
  name: string;
  pubchem_cid?: string | null;
  broad_id?: string | null;
  chembl_id?: string | null;
  smiles?: string | null;
  inchikey?: string | null;
  molecular_formula?: string | null;
  molecular_weight?: number | null;
  iupac_name?: string | null;
}

/** Renders "C29H31N7O" with the digits subscripted. */
function Formula({ formula }: { formula: string }) {
  return (
    <span className="font-mono">
      {formula.split(/(\d+)/).map((part, i) =>
        /^\d+$/.test(part)
          ? <sub key={i} className="text-[0.7em]">{part}</sub>
          : <span key={i}>{part}</span>,
      )}
    </span>
  );
}

function Field({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div className="flex flex-col gap-0.5 min-w-0">
      <span className="text-[10px] uppercase tracking-widest text-muted-foreground font-semibold">
        {label}
      </span>
      <span className="text-sm text-foreground break-words">{children}</span>
    </div>
  );
}

const DASH = <span className="text-muted-foreground/60">—</span>;

/**
 * Summary of the currently selected drug: structure on the left, identity on the right.
 *
 * The 3D view is the drug's ACTUAL bound conformer, sliced out of the co-folded AF3 complex
 * (ligandOnly hides the protein) — not a generic PubChem conformer. The 2D skeleton is drawn
 * locally from SMILES. Neither makes an external request.
 */
export default function DrugSummaryPanel({
  drug,
  complexCifUrl,
}: {
  drug?: DrugSummary | null;
  complexCifUrl?: string | null;
}) {
  if (!drug) return null;

  const pubchemUrl = drug.pubchem_cid
    ? `https://pubchem.ncbi.nlm.nih.gov/compound/${drug.pubchem_cid}`
    : null;
  // InChIKey is the stable cross-database handle — resolve it rather than guessing a
  // ChEMBL id we do not have (chembl_id is still unpopulated).
  const chemblUrl = drug.chembl_id
    ? `https://www.ebi.ac.uk/chembl/compound_report_card/${drug.chembl_id}/`
    : drug.inchikey
      ? `https://www.ebi.ac.uk/unichem/compoundsources?type=inchikey&compound=${drug.inchikey}`
      : null;

  return (
    <Card className="overflow-hidden">
      <CardHeader className="py-4">
        <CardTitle className="text-lg font-semibold flex items-center gap-2 flex-wrap">
          <FlaskConical className="w-4 h-4 text-emerald-500 shrink-0" />
          <span className="truncate">{drug.name}</span>
          {drug.molecular_formula && (
            <Badge variant="secondary" className="text-[10px] font-normal">
              <Formula formula={drug.molecular_formula} />
            </Badge>
          )}
        </CardTitle>
      </CardHeader>

      <CardContent className="pt-0">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          {/* LEFT — 2D skeleton above the bound 3D conformer */}
          <div className="flex flex-col gap-2">
            <div className="rounded-lg border border-border/40 bg-muted/20 h-[150px] overflow-hidden">
              <MoleculeSkeleton smiles={drug.smiles} className="w-full h-full" />
            </div>

            <div className="rounded-lg border border-border/40 bg-muted/20 h-[190px] overflow-hidden relative">
              {complexCifUrl ? (
                <StructureViewer
                  key={`lig-${drug.id}-${complexCifUrl}`}
                  cifUrl={complexCifUrl}
                  ligandBallAndStick
                  ligandOnly
                  primaryLabel={drug.name}
                />
              ) : (
                <div className="h-full flex items-center justify-center text-xs text-muted-foreground">
                  No co-folded structure for this drug
                </div>
              )}
            </div>
            <p className="text-[10px] text-muted-foreground leading-snug">
              Top: 2D skeletal formula drawn from SMILES. Bottom: the drug&apos;s bound conformer,
              taken from its co-folded AlphaFold 3 complex (protein hidden).
            </p>
          </div>

          {/* RIGHT — identity & chemistry */}
          <div className="flex flex-col gap-3">
            <div className="grid grid-cols-2 gap-3">
              <Field label="Molecular Formula">
                {drug.molecular_formula ? <Formula formula={drug.molecular_formula} /> : DASH}
              </Field>
              <Field label="Molar Mass">
                {drug.molecular_weight != null
                  ? <span className="font-mono">{drug.molecular_weight.toFixed(2)} g/mol</span>
                  : DASH}
              </Field>
            </div>

            <div className="h-px w-full bg-border/50" />

            <div className="grid grid-cols-2 gap-3">
              <Field label="PubChem CID">
                {pubchemUrl ? (
                  <a href={pubchemUrl} target="_blank" rel="noopener noreferrer"
                     className="font-mono text-primary/90 hover:text-primary hover:underline inline-flex items-center gap-1">
                    {drug.pubchem_cid}
                    <ExternalLink className="w-2.5 h-2.5 opacity-70" />
                  </a>
                ) : DASH}
              </Field>
              <Field label={drug.chembl_id ? "ChEMBL" : "UniChem"}>
                {chemblUrl ? (
                  <a href={chemblUrl} target="_blank" rel="noopener noreferrer"
                     className="font-mono text-primary/90 hover:text-primary hover:underline inline-flex items-center gap-1">
                    {drug.chembl_id ?? "lookup"}
                    <ExternalLink className="w-2.5 h-2.5 opacity-70" />
                  </a>
                ) : DASH}
              </Field>
            </div>

            <Field label="InChIKey">
              {drug.inchikey
                ? <span className="font-mono text-xs break-all">{drug.inchikey}</span>
                : DASH}
            </Field>

            {drug.broad_id && (
              <Field label="Broad ID">
                <span className="font-mono text-xs">{drug.broad_id}</span>
              </Field>
            )}

            <Field label="SMILES">
              {drug.smiles
                ? <span className="font-mono text-[11px] break-all leading-relaxed block max-h-14 overflow-y-auto pr-1">
                    {drug.smiles}
                  </span>
                : DASH}
            </Field>

            {drug.iupac_name && (
              <Field label="IUPAC Name">
                <span className="text-xs text-muted-foreground leading-relaxed block max-h-14 overflow-y-auto pr-1">
                  {drug.iupac_name}
                </span>
              </Field>
            )}
          </div>
        </div>
      </CardContent>
    </Card>
  );
}
