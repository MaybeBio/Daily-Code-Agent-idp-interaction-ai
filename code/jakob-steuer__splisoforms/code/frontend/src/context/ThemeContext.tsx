"use client";

import React, { createContext, useContext, useEffect, useState } from "react";

type Theme = "light" | "dark";

interface ThemeContextType {
  theme: Theme;
  toggleTheme: () => void;
  setTheme: (theme: Theme) => void;
}

const ThemeContext = createContext<ThemeContextType | undefined>(undefined);

/** Light is the default — explicit, both for first-time visitors and for any
 *  cleared-localStorage state. The actual applied class is set by the
 *  blocking inline script in `<head>` (see `app/layout.tsx`), and we mirror
 *  it here so the very first React render carries the correct value. */
const DEFAULT_THEME: Theme = "light";

function readInitialTheme(): Theme {
  // SSR: must be deterministic, can't read the DOM. Match the inline-script
  // default so the SSR markup lines up with what the script will produce.
  if (typeof document === "undefined") return DEFAULT_THEME;
  return document.documentElement.classList.contains("dark") ? "dark" : "light";
}

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  // Lazy initializer reads the class that the blocking <head> script has
  // already set, so theme === actual visible state from the very first render.
  // No more "useState('dark') -> useEffect resets to 'light'" bounce.
  const [theme, setThemeState] = useState<Theme>(readInitialTheme);

  const setTheme = (newTheme: Theme) => {
    setThemeState(newTheme);
    try {
      localStorage.setItem("spliso-theme", newTheme);
    } catch {
      /* private mode / disabled storage — ignore */
    }
    const root = document.documentElement;
    root.classList.toggle("dark", newTheme === "dark");
    root.style.colorScheme = newTheme;
  };

  const toggleTheme = () => {
    setTheme(theme === "dark" ? "light" : "dark");
  };

  return (
    <ThemeContext.Provider value={{ theme, toggleTheme, setTheme }}>
      {children}
    </ThemeContext.Provider>
  );
}

export function useTheme() {
  const context = useContext(ThemeContext);
  if (context === undefined) {
    throw new Error("useTheme must be used within a ThemeProvider");
  }
  return context;
}
