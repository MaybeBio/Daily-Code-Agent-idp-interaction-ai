"use client";

import { usePathname } from "next/navigation";
import Link from "next/link";
import { cn } from "@/lib/utils";
import {
  Book,
  Layers,
  Dna,
  Activity,
  ShieldAlert,
  Zap,
  Compass,
  Scissors,
  Network,
  ArrowLeft,
  ChevronRight,
  Menu,
  X,
  Pill
} from "lucide-react";
import { Button } from "@/components/ui/button";
import { useState } from "react";

const navigation = [
  { name: "Overview", href: "/documentation", icon: Book },
  { name: "Protein Structure", href: "/documentation/protein-structure", icon: Layers },
  { name: "Functional Domains", href: "/documentation/functional-domains", icon: Dna },
  { name: "NMD & Transcript Fate", href: "/documentation/nmd-transcript-fate", icon: ShieldAlert },
  { name: "Neoantigens", href: "/documentation/neoantigens", icon: Zap },
  { name: "PTM Accessibility", href: "/documentation/ptm-accessibility", icon: Activity },
  { name: "Interactions & Expression", href: "/documentation/interactions-expression", icon: Network },
  { name: "Drug Interactions", href: "/documentation/drug-interactions", icon: Pill },
  { name: "Discovery & Browsing", href: "/documentation/discovery", icon: Compass },
];

export default function DocumentationLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  const pathname = usePathname();
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  return (
    <div className="flex min-h-screen bg-background pt-16">
      {/* Mobile Sidebar Toggle */}
      <div className="lg:hidden fixed bottom-6 right-6 z-50">
        <Button 
          size="icon" 
          onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
          className="rounded-full h-12 w-12 shadow-xl border border-primary/20"
        >
          {mobileMenuOpen ? <X className="w-6 h-6" /> : <Menu className="w-6 h-6" />}
        </Button>
      </div>

      {/* Sidebar Overlay */}
      {mobileMenuOpen && (
        <div 
          className="fixed inset-0 bg-background/80 backdrop-blur-sm z-40 lg:hidden"
          onClick={() => setMobileMenuOpen(false)}
        />
      )}

      {/* Navigation Sidebar */}
      <aside className={cn(
        "fixed lg:sticky top-16 h-[calc(100vh-4rem)] w-64 border-r border-border/40 bg-muted/30 backdrop-blur-md z-40 transition-transform lg:translate-x-0 overflow-y-auto",
        mobileMenuOpen ? "translate-x-0" : "-translate-x-full"
      )}>
        <div className="p-6 space-y-6">
          <div className="space-y-1">
            <Link href="/" className="flex items-center gap-2 text-muted-foreground hover:text-primary transition-colors text-xs font-bold uppercase tracking-widest mb-6 px-3">
              <ArrowLeft className="w-3 h-3" />
              Main Dashboard
            </Link>
          </div>

          <nav className="space-y-1.5">
            {navigation.map((item) => {
              const active = pathname === item.href;
              return (
                <Link
                  key={item.name}
                  href={item.href}
                  onClick={() => setMobileMenuOpen(false)}
                  className={cn(
                    "group flex items-center gap-3 px-3 py-2.5 text-sm font-medium rounded-lg transition-all duration-200 border border-transparent",
                    active 
                      ? "bg-primary/10 text-primary border-primary/20 shadow-sm" 
                      : "text-muted-foreground hover:bg-muted hover:text-foreground hover:border-border/60"
                  )}
                >
                  <item.icon className={cn(
                    "w-4 h-4 transition-colors",
                    active ? "text-primary" : "text-muted-foreground group-hover:text-foreground"
                  )} />
                  {item.name}
                  {active && <ChevronRight className="ml-auto w-3 h-3 text-primary/60" />}
                </Link>
              );
            })}
          </nav>
        </div>
      </aside>

      {/* Main Content Area.
          Content width adapts to the viewport in bounded steps rather than a single
          fixed cap: comfortable prose measure on laptops, more room for the card
          grids / tables on large monitors, without going full-bleed (which would
          stretch text lines past a readable length). Dial the 2xl step down to
          max-w-6xl if the prose feels too wide on very large screens. */}
      <main className="flex-1 min-w-0 flex flex-col items-center">
        <div className="w-full max-w-4xl lg:max-w-5xl xl:max-w-6xl 2xl:max-w-7xl px-6 lg:px-12 py-12 lg:py-16">
          {children}
        </div>
      </main>
    </div>
  );
}
