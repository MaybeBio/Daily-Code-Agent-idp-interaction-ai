import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Emit source maps for the production client bundle so runtime errors
  // surface real file:line locations in the browser console instead of
  // post-minification gibberish like "1d410b06af187467.js:5:95802".
  // Trade-off: ~5-10% bigger output, but invaluable for diagnosing the
  // class of crashes that only manifest in `next build` mode.
  productionBrowserSourceMaps: true,

  // Mol* (molstar) v5.x compiles its TypeScript namespaces to the legacy
  //   var X; (function (X) { X.Provider = ...; })(X || (X = {}));
  //   export { X };
  // pattern, then mmcif.js circularly imports modules like anisotropic.js
  // and reads `X.Provider` at top-level (line 54 et al). Turbopack's
  // production output mishandles the live ESM binding through that
  // circular cycle and resolves `X` to undefined → "Cannot read properties
  // of undefined (reading 'Provider')". We avoid the issue entirely by
  // staying on the Next.js 15 webpack production pipeline, which handles
  // circular ESM correctly. Keeping molstar in transpilePackages so any
  // future Turbopack experiment will at least re-process it through SWC.
  transpilePackages: ["molstar"],
  eslint: {
    ignoreDuringBuilds: true,
  },
};

export default nextConfig;

// bump: full container rebuild after podman system reset, 2026-06-12
