-- ============================================================================
-- Read-only Postgres role for the public REST API (/api/v1).
--
-- Run once as a superuser / DB owner on the splisoforms database:
--     psql -d splisoforms -f backend/scripts/setup_readonly_role.sql
--
-- Then point the backend at it:
--     export PUBLIC_DATABASE_URL=postgresql://splisoforms_public:CHANGE_ME@localhost:5432/splisoforms
--
-- Defense in depth: even if a bug or injection reaches the public surface, this
-- role can only SELECT, and a statement_timeout caps run-away queries.
-- ============================================================================

-- 1. The role (pick a strong password and mirror it in PUBLIC_DATABASE_URL).
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = 'splisoforms_public') THEN
        CREATE ROLE splisoforms_public LOGIN PASSWORD 'CHANGE_ME';
    END IF;
END
$$;

-- 2. No write ability anywhere; explicit SELECT on existing + future tables.
REVOKE ALL ON ALL TABLES    IN SCHEMA public FROM splisoforms_public;
REVOKE ALL ON ALL SEQUENCES IN SCHEMA public FROM splisoforms_public;

GRANT CONNECT ON DATABASE splisoforms TO splisoforms_public;
GRANT USAGE  ON SCHEMA public         TO splisoforms_public;
GRANT SELECT ON ALL TABLES IN SCHEMA public TO splisoforms_public;

-- Future tables created by migrations also become readable automatically.
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO splisoforms_public;

-- 3. Cap run-away queries (10s) for this role specifically.
ALTER ROLE splisoforms_public SET statement_timeout = '10s';

-- 4. (Optional, stronger) Revoke read on internal/IP-sensitive tables so the
--    public role physically cannot read them even if an endpoint slips through.
REVOKE SELECT ON prediction_jobs FROM splisoforms_public;
REVOKE SELECT ON api_keys        FROM splisoforms_public;
