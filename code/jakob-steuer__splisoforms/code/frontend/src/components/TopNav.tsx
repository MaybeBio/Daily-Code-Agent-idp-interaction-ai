"use client";

import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Dna } from "lucide-react";
import { Button, buttonVariants } from "@/components/ui/button";
import Link from "next/link";
import { SearchCatalog } from "@/components/SearchCatalog";
import { cn } from "@/lib/utils";

import { useState, useEffect } from "react";
import { usePathname } from "next/navigation";
import { ThemeToggle } from "@/components/ThemeToggle";

// Primary "lens" destinations — each a proteome-wide angle that drills into
// the shared structure viewer.
const lenses = [
  { name: "Structure",       href: "/structure" },
  { name: "Splicing Impact", href: "/explore/impact" },
  { name: "Neoepitopes",     href: "/explore/neoepitopes" },
  { name: "Drugs",           href: "/explore/drugs" },
];

export function TopNav() {
  const [mounted, setMounted] = useState(false);
  const pathname = usePathname();

  useEffect(() => {
    setMounted(true);
  }, []);

  return (
    <nav className="relative z-[100] w-full max-w-[1800px] mx-auto px-6 py-4 flex flex-col md:flex-row justify-between items-center border-b border-border/50 gap-4 lg:gap-8">
      <div className="shrink-0 flex items-center">
        <Link href="/" className={cn("flex items-center gap-3", mounted ? "cursor-pointer group" : "")}>
          <Dna className={cn("w-8 h-8 text-primary", mounted ? "animate-pulse group-hover:scale-110 transition-transform" : "opacity-50")} />
          <span className={cn("font-bold text-xl tracking-widest", mounted ? "text-transparent bg-clip-text bg-gradient-to-r from-primary to-secondary" : "text-slate-500")}>
            SPLISOFORMS
          </span>
        </Link>
      </div>

      <div className="flex-1 min-w-0 w-full max-w-2xl order-3 md:order-none mt-2 md:mt-0">
         <SearchCatalog compact={true} />
      </div>


      <div className="flex gap-1.5 shrink-0 items-center overflow-x-auto pb-2 md:pb-0 w-full md:w-auto justify-center md:justify-end">
        {/* Primary lenses */}
        {lenses.map((lens) => {
          const active = pathname === lens.href || pathname.startsWith(lens.href + "/");
          return (
            <Link
              key={lens.href}
              href={lens.href}
              className={cn(
                buttonVariants({ variant: active ? "default" : "ghost" }),
                "whitespace-nowrap",
                !active && "hover:text-primary"
              )}
            >
              {lens.name}
            </Link>
          );
        })}

        {/* Divider */}
        <span className="hidden md:inline-block w-px h-5 bg-border/60 mx-1" />

        {/* Utility */}
        <Link href="/documentation" className={cn(buttonVariants({ variant: "ghost", size: "sm" }), "text-muted-foreground hover:text-primary")}>
          Docs
        </Link>
        <Link href="/api-docs" className={cn(buttonVariants({ variant: "ghost", size: "sm" }), "text-muted-foreground hover:text-secondary")}>
          API
        </Link>
        <ThemeToggle />
      </div>
    </nav>
  );
}
