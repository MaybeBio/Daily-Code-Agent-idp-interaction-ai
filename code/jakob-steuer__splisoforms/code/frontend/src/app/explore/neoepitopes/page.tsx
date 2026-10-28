"use client";

import { Suspense } from "react";
import { useRouter } from "next/navigation";
import { ArrowLeft } from "lucide-react";
import { Button } from "@/components/ui/button";
import GlobalNeoantigenBrowser from "@/components/GlobalNeoantigenBrowser";

function NeoepitopesInner() {
  const router = useRouter();

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-[1700px] mx-auto p-6 md:p-10">
        <header className="flex flex-col md:flex-row justify-between items-start md:items-center gap-4 mb-8">
          <div className="flex gap-4 items-center">
            <Button variant="ghost" size="icon" onClick={() => router.push("/")} className="hover:bg-rose-500/10">
              <ArrowLeft className="w-5 h-5" />
            </Button>
            <div>
              <h1 className="text-4xl font-bold tracking-tight bg-clip-text text-transparent bg-gradient-to-r from-rose-400 via-foreground/80 to-secondary">
                Neoepitope Atlas
              </h1>
              <p className="text-sm text-muted-foreground mt-1">
                Globally ranked splicing-derived neoantigens across all isoforms, scored by an
                MHCflurry · NetMHCpan · BigMHC ensemble.
              </p>
            </div>
          </div>

        </header>

        <GlobalNeoantigenBrowser />
      </div>
    </div>
  );
}

export default function NeoepitopesPage() {
  return (
    <Suspense fallback={<div className="text-center py-24 animate-pulse text-rose-400 tracking-widest uppercase text-xs">Initializing Neoepitope Atlas...</div>}>
      <NeoepitopesInner />
    </Suspense>
  );
}
