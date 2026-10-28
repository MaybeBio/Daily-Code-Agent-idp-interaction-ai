import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";
import { TopNav } from "@/components/TopNav";
import { Footer } from "@/components/Footer";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "SPLISOFORMS — Splicing Proteomics Knowledgebase",
  description: "High-performance exploration of isoform/structure data.",
};

import { ThemeProvider } from "@/context/ThemeContext";

// Blocking script — runs synchronously in <head> before React hydrates and
// before any pixel paints. Reads the persisted theme from localStorage and
// applies the dark class + colorScheme to <html>. Defaults to "light" when
// nothing is stored. Without this the page would flash the wrong theme on
// every reload (and Mol* would init its canvas with the wrong background).
const themeBootstrap = `
(function() {
  try {
    var t = localStorage.getItem('spliso-theme');
    var theme = (t === 'dark' || t === 'light') ? t : 'light';
    var root = document.documentElement;
    if (theme === 'dark') {
      root.classList.add('dark');
      root.style.colorScheme = 'dark';
    } else {
      root.classList.remove('dark');
      root.style.colorScheme = 'light';
    }
  } catch (e) {}
})();
`;

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <script dangerouslySetInnerHTML={{ __html: themeBootstrap }} />
        {/* Preload both home-hero backgrounds so the late image-load doesn't
            trigger a repaint of the hero overlay on first interactive frame
            (was visible as a flicker of the top search-bar and release pill
            on the production server, not in dev where the cache was warm). */}
        <link rel="preload" as="image" href="/background-light.webp" />
        <link rel="preload" as="image" href="/background-dark.webp" />
      </head>
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased min-h-screen flex flex-col`}
        suppressHydrationWarning
      >
        <ThemeProvider>
          <TopNav />
          <main className="flex-1 flex flex-col">
            {children}
          </main>
          <Footer />
        </ThemeProvider>
      </body>
    </html>
  );
}
