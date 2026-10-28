/**
 * Pick the correct API base URL for the current execution context.
 *
 * Why this exists
 * ---------------
 * Next.js Server Components run `fetch(...)` INSIDE the container at SSR
 * time. If that fetch targets the same public domain the browser uses
 * (`https://splisoforms.org`), the container has to resolve its own
 * public IP, hairpin-NAT back through the host firewall, and validate
 * the TLS cert from inside its own user-namespace. On the FHNW server
 * this either times out or fails CA validation — the SSR page render
 * throws, returns 500, and the client-side hydration crashes with the
 * confusing `Cannot read properties of undefined (reading 'Provider')`.
 *
 * The fix: split SSR and CSR URLs.
 *
 *   * SSR   → http://127.0.0.1:8000      (server-internal; no TLS, no DNS)
 *   * CSR   → "" (empty, => relative URLs) — same-origin through nginx
 *
 * Local dev override:
 *   * SSR   → NEXT_PUBLIC_API_URL if set (`http://localhost:8000` for the
 *             dev backend; not bound to a container)
 *   * CSR   → NEXT_PUBLIC_API_URL if set, else "" (relative)
 *
 * Returns a base URL that callers concatenate with absolute API paths:
 *
 *     fetch(`${getApiUrl()}/api/genes/${id}`)
 *
 * If the base is empty (CSR + same-origin), the resulting fetch URL is
 * `/api/genes/...` — the browser resolves it against the current origin.
 */
export function getApiUrl(): string {
  const isServer = typeof window === "undefined";

  if (isServer) {
    // Container-side fetch. Prefer the explicit override so we can point
    // a staging build at a remote dev backend if we ever need to.
    return process.env.INTERNAL_API_URL || "http://127.0.0.1:8000";
  }

  // Browser-side fetch. NEXT_PUBLIC_API_URL is set at build time
  // (Dockerfile pins it for prod, .env.local pins it for local dev).
  // If the value is the production public domain, an empty string is
  // safer (same-origin via nginx — avoids CORS concerns).
  const pub = process.env.NEXT_PUBLIC_API_URL || "";
  if (pub.startsWith("https://splisoforms.org") || pub === "/") {
    // Same-origin: let the browser resolve against the current URL.
    return "";
  }
  return pub;
}
