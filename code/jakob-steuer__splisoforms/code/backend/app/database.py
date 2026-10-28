import os
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker

SQLALCHEMY_DATABASE_URL = os.getenv("DATABASE_URL", "postgresql://splisoforms_user:splisoforms_dev@localhost:5432/splisoforms")

engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    pool_size=20,
    max_overflow=10,
    pool_pre_ping=True,
    pool_recycle=3600,
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()


# ── Public API read-only session ───────────────────────────────────────────────
# The /api/v1 public endpoints connect through a separate, SELECT-only Postgres
# role (see scripts/setup_readonly_role.sql) with a statement_timeout, so a bug
# or injection in the public surface can neither write nor run away with a
# long-running query.  Falls back to the main URL (with a warning) when
# PUBLIC_DATABASE_URL is unset, so dev still works without the role provisioned.
PUBLIC_DATABASE_URL = os.getenv("PUBLIC_DATABASE_URL")
if not PUBLIC_DATABASE_URL:
    import warnings
    warnings.warn(
        "PUBLIC_DATABASE_URL not set — the public API will use the main "
        "(write-capable) DB role. Provision the read-only role for production.",
        RuntimeWarning,
    )
    PUBLIC_DATABASE_URL = SQLALCHEMY_DATABASE_URL

public_engine = create_engine(
    PUBLIC_DATABASE_URL,
    pool_size=10,
    max_overflow=5,
    pool_pre_ping=True,
    pool_recycle=3600,
)
PublicSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=public_engine)


def create_worker_engine(pool_size: int = 2, max_overflow: int = 1):
    """Create a lightweight engine for use in ProcessPoolExecutor workers.

    Each forked process must create its own engine — sharing the parent's
    engine across process boundaries corrupts connections. Use small pool
    sizes to avoid exhausting the database's connection limit.
    """
    return create_engine(
        SQLALCHEMY_DATABASE_URL,
        pool_size=pool_size,
        max_overflow=max_overflow,
        pool_pre_ping=True,
        pool_recycle=3600,
    )


def get_db():
    db = SessionLocal()
    try:
        yield db
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def get_public_db():
    """Read-only session for the public /api/v1 endpoints (SELECT-only role)."""
    db = PublicSessionLocal()
    try:
        yield db
    finally:
        db.close()
