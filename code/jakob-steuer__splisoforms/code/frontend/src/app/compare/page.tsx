"use client";

import { Suspense } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import { ArrowLeft, Loader2 } from "lucide-react";
import { Button } from "@/components/ui/button";
import CompareViewClient from "@/components/CompareViewClient";
import CompareBuilder from "@/components/CompareBuilder";

function ComparisonWrapper() {
  const searchParams = useSearchParams();
  const idsParam = searchParams.get("ids");
  const router = useRouter();

  if (!idsParam) {
    return <CompareBuilder />;
  }

  const ids = Array.from(new Set(
    idsParam.split(",").map(id => parseInt(id, 10)).filter(id => !isNaN(id))
  )).slice(0, 5);

  if (ids.length < 2) {
    return (
      <div className="flex flex-col items-center justify-center min-h-[60vh] text-center p-10">
        <h2 className="text-2xl font-bold text-foreground mb-2">Need at least two isoforms</h2>
        <p className="text-muted-foreground mb-6"> The comparison view needs between 2 and 5 isoforms. </p>
        <Button onClick={() => router.back()}>Go Back</Button>
      </div>
    );
  }

  return <CompareViewClient ids={ids} />;
}

export default function ComparePage() {
  const router = useRouter();

  return (
    <div className="min-h-screen bg-background text-foreground pb-20">
      <header className="border-b border-primary/20 bg-card/50 backdrop-blur sticky top-0 z-50">
        <div className="w-full max-w-[1800px] mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-4">
            <Button variant="ghost" size="icon" onClick={() => router.back()} className="hover:bg-muted text-foreground">
              <ArrowLeft className="w-5 h-5" />
            </Button>
            <h1 className="text-xl font-bold bg-gradient-to-r from-primary to-secondary bg-clip-text text-transparent">
              Isoform Comparison Mode
            </h1>
          </div>
          <div className="text-xs text-muted-foreground uppercase tracking-widest font-medium">
            Side-by-Side Analysis
          </div>
        </div>
      </header>

      <main className="w-full max-w-[1800px] mx-auto px-6 py-8">
        <Suspense fallback={
          <div className="flex flex-col items-center justify-center min-h-[50vh] gap-4">
            <Loader2 className="w-10 h-10 text-primary animate-spin" />
            <p className="text-muted-foreground">Loading isoforms for comparison...</p>
          </div>
        }>
          <ComparisonWrapper />
        </Suspense>
      </main>
    </div>
  );
}
