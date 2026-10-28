import { notFound } from "next/navigation";
import GeneHeader from "./components/GeneHeader";
import GeneDashboardClient from "./components/GeneDashboardClient";
import GeneErrorFallback from "./components/GeneErrorFallback";
import { getApiUrl } from "@/lib/api";

export default async function GeneOverviewPage({ params }: { params: Promise<{ id: string }> }) {
  const { id } = await params;
  const apiUrl = getApiUrl();
  
  let gene;
  try {
    const res = await fetch(`${apiUrl}/api/genes/${id}`, { 
      cache: "no-store",
      next: { revalidate: 0 } 
    });
    
    if (!res.ok) {
      if (res.status === 404) notFound();
      throw new Error(`Failed to fetch gene data: ${res.status}`);
    }
    
    gene = await res.json();
  } catch (error: any) {
    console.error("Gene fetch error:", error);
    return <GeneErrorFallback id={id} error={error.message || "Unknown error"} />;
  }

  return (
    <div className="container mx-auto py-8 px-4 flex flex-col gap-8 max-w-[1600px] animate-in fade-in duration-500">
      <GeneHeader gene={gene} />
      <GeneDashboardClient gene={gene} />
    </div>
  );
}
