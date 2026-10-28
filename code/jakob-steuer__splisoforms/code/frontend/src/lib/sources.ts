// Long DB source names (the `sources` table) → concise display labels.
// The full name stays available as a tooltip; the /documentation page carries the
// detailed dataset description. Keep this map in sync with the `sources` table.
const SHORT_LABELS: Record<string, string> = {
  "ccRCC PacBio MAS-seq (Karakulak 2025)": "ccRCC (Karakulak 2025)",
  "Breast Cancer Nanopore ESPRESSO (Head 2026)": "Breast Cancer (Head 2026)",
};

/** Concise label for a dataset source name; falls back to the name itself. */
export function sourceShortLabel(name?: string | null): string {
  if (!name) return "";
  return SHORT_LABELS[name] ?? name;
}
