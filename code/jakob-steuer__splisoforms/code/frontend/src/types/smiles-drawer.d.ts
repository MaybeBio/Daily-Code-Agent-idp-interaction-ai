/**
 * Minimal typings for `smiles-drawer` (v2), which ships no .d.ts of its own.
 * Only the surface we use: parse a SMILES string, then draw it into an <svg>.
 */
declare module "smiles-drawer" {
  interface DrawerOptions {
    width?: number;
    height?: number;
    padding?: number;
    bondThickness?: number;
    bondLength?: number;
    fontSizeLarge?: number;
    fontSizeSmall?: number;
    terminalCarbons?: boolean;
    explicitHydrogens?: boolean;
    themes?: Record<string, Record<string, string>>;
  }

  class SvgDrawer {
    constructor(options?: DrawerOptions);
    /** `target` may be an SVGElement or an element id. `theme` keys into options.themes. */
    draw(tree: unknown, target: SVGElement | string, theme?: string, infoOnly?: boolean): void;
  }

  const SmilesDrawer: {
    SvgDrawer: typeof SvgDrawer;
    /** Async parse; calls `onError` for invalid SMILES rather than throwing. */
    parse(
      smiles: string,
      onSuccess: (tree: unknown) => void,
      onError?: (err: unknown) => void,
    ): void;
  };

  export default SmilesDrawer;
}
