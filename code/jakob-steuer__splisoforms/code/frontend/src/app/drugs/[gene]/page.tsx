import { notFound } from "next/navigation";
import { getApiUrl } from "@/lib/api";
import DrugGeneClient from "./DrugGeneClient";

export default async function DrugGenePage({ params }: { params: Promise<{ gene: string }> }) {
  const { gene } = await params;
  const apiUrl = getApiUrl();

  let geneData: any;
  let drugsData: any;
  try {
    const [gRes, dRes] = await Promise.all([
      fetch(`${apiUrl}/api/genes/${gene}`, { cache: "no-store" }),
      fetch(`${apiUrl}/api/gene/${gene}/drugs`, { cache: "no-store" }),
    ]);
    if (gRes.status === 404 || dRes.status === 404) notFound();
    if (!gRes.ok || !dRes.ok) throw new Error(`Failed to load drug data (${gRes.status}/${dRes.status})`);
    geneData = await gRes.json();
    drugsData = await dRes.json();
  } catch (err: any) {
    return (
      <div className="container mx-auto py-16 px-4 text-center text-muted-foreground">
        Could not load drug data for {gene}: {err.message}
      </div>
    );
  }

  if (!drugsData.complexes?.length) {
    return (
      <div className="container mx-auto py-16 px-4 text-center text-muted-foreground">
        No co-folded drug complexes for <span className="font-bold text-foreground">{geneData.symbol}</span>.
      </div>
    );
  }

  return <DrugGeneClient gene={geneData} drugs={drugsData} />;
}
