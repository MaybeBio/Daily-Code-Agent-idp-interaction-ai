import os
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from app.api.bio import router as bio_router
from app.api.v1 import router as v1_router
from app.database import engine, Base

# Create DB tables (if not already existing)
# In production, use Alembic migrations instead
# Base.metadata.create_all(bind=engine)

app = FastAPI(title="SPLISOFORMS API", version="1.1.0")

# Enable GZip compression (crucial for large PAE JSON matrices)
app.add_middleware(GZipMiddleware, minimum_size=1000)

# Configure CORS.
#
# Origins are read from the SPLISOFORMS_CORS_ORIGINS env var (comma-separated
# list, e.g. "https://your-domain.org,https://www.your-domain.org").
# Falls back to a dev-only default when the var is unset so local development
# keeps working without ceremony. In production set this in the server's .env
# to the actual public origins — never leave it on the wildcard once the API
# is reachable from the internet.
#
# allow_credentials must be False when using a wildcard origin (browsers
# reject the combination per CORS spec). We keep credentials True in the
# explicit-origin path and flip to False if a wildcard slips through.
_cors_env = os.getenv("SPLISOFORMS_CORS_ORIGINS", "").strip()
if _cors_env:
    _cors_origins = [o.strip() for o in _cors_env.split(",") if o.strip()]
else:
    # Dev fallback — three localhost ports covering the typical Next.js +
    # FastAPI dev loop. Production deployments MUST set the env var.
    _cors_origins = [
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:8000",
    ]

_uses_wildcard = "*" in _cors_origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=not _uses_wildcard,
    allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
    allow_headers=["Authorization", "Content-Type", "Accept", "X-Requested-With"],
)

# Mount static folder for structural files
static_dir = os.path.join(os.path.dirname(__file__), "static")
if not os.path.exists(static_dir):
    os.makedirs(os.path.join(static_dir, "data"), exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────
# Transparent gzip fallback for AF3 PAE / confidences JSON files.
#
# After the pre-server-transfer cleanup pass, all `*_confidences.json` files
# are stored on disk as `*_confidences.json.gz`.  The browser still requests
# the original `.json` URL — we serve the `.gz` body with the right
# Content-Encoding header so it decompresses transparently in fetch().
#
# Must be registered BEFORE `app.mount("/static", ...)` because FastAPI/
# Starlette routes are matched in declaration order.
# ──────────────────────────────────────────────────────────────────────────

@app.get("/static/{path:path}_confidences.json")
async def serve_confidences_json(path: str):
    """Serve <file>_confidences.json, falling back to <file>_confidences.json.gz."""
    base = os.path.join(static_dir, f"{path}_confidences.json")
    if os.path.exists(base):
        return FileResponse(base, media_type="application/json")
    gz = base + ".gz"
    if os.path.exists(gz):
        # Returning the gzipped bytes with the Content-Encoding header lets
        # the browser/fetch decompress on the wire — same as nginx
        # gzip_static.  Mol* / our PAEMatrix fetcher work unchanged.
        return FileResponse(
            gz,
            media_type="application/json",
            headers={"Content-Encoding": "gzip"},
        )
    raise HTTPException(status_code=404, detail="confidences file not found")


app.mount("/static", StaticFiles(directory=static_dir, follow_symlink=True), name="static")

# Include the bio API endpoints (internal — used by the frontend)
app.include_router(bio_router, prefix="/api", tags=["bio"])

# Public, curated, read-only API (API-key auth + per-key rate limit)
app.include_router(v1_router, prefix="/api/v1", tags=["public-v1"])

@app.get("/health")
async def health_check():
    return {"status": "ok", "message": "SPLISOFORMS API is running"}

# bump: full container rebuild after podman system reset, 2026-06-12
