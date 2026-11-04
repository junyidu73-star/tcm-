-- PostgreSQL production schema for the evidence-first workflow.
-- SQLite is used for this single-computer demonstration run; this DDL is the
-- deployable relational schema and intentionally does not require a graph DB.
CREATE TABLE IF NOT EXISTS entity (
  entity_id TEXT PRIMARY KEY,
  entity_type TEXT NOT NULL,
  namespace TEXT,
  canonical_name TEXT NOT NULL,
  species TEXT,
  created_at TIMESTAMPTZ NOT NULL DEFAULT now(),
  metadata JSONB NOT NULL DEFAULT '{}'::jsonb
);
CREATE TABLE IF NOT EXISTS entity_alias (
  alias_id BIGSERIAL PRIMARY KEY,
  entity_id TEXT NOT NULL REFERENCES entity(entity_id),
  source TEXT NOT NULL,
  source_id TEXT,
  raw_name TEXT NOT NULL,
  normalized_name TEXT,
  mapping_status TEXT NOT NULL CHECK (mapping_status IN ('exact','synonym','ambiguous','unmatched','manual')),
  mapping_algorithm TEXT,
  mapping_version TEXT,
  UNIQUE(entity_id, source, raw_name)
);
CREATE TABLE IF NOT EXISTS study (
  study_id TEXT PRIMARY KEY,
  source_database TEXT NOT NULL,
  source_accession TEXT NOT NULL,
  title TEXT,
  assay TEXT,
  data_modality TEXT,
  species TEXT,
  tissue TEXT,
  cell_type TEXT,
  experimental_design JSONB NOT NULL DEFAULT '{}'::jsonb,
  source_version TEXT,
  source_url TEXT,
  retrieved_at TIMESTAMPTZ NOT NULL DEFAULT now(),
  UNIQUE(source_database, source_accession)
);
CREATE TABLE IF NOT EXISTS sample (
  sample_id TEXT PRIMARY KEY,
  study_id TEXT NOT NULL REFERENCES study(study_id),
  group_label TEXT,
  diagnosis TEXT,
  species TEXT,
  tissue TEXT,
  cell_type TEXT,
  intervention_id TEXT,
  dose NUMERIC,
  dose_unit TEXT,
  duration NUMERIC,
  duration_unit TEXT,
  batch TEXT,
  metadata JSONB NOT NULL DEFAULT '{}'::jsonb
);
CREATE TABLE IF NOT EXISTS intervention (
  intervention_id TEXT PRIMARY KEY,
  intervention_type TEXT NOT NULL,
  entity_id TEXT REFERENCES entity(entity_id),
  raw_name TEXT NOT NULL,
  dose NUMERIC,
  dose_unit TEXT,
  route TEXT,
  duration NUMERIC,
  duration_unit TEXT,
  preparation TEXT,
  metadata JSONB NOT NULL DEFAULT '{}'::jsonb
);
CREATE TABLE IF NOT EXISTS evidence_edge (
  evidence_id TEXT PRIMARY KEY,
  subject_id TEXT NOT NULL REFERENCES entity(entity_id),
  predicate TEXT NOT NULL,
  object_id TEXT NOT NULL REFERENCES entity(entity_id),
  source_database TEXT NOT NULL,
  source_id TEXT NOT NULL,
  evidence_type TEXT NOT NULL,
  evidence_level TEXT,
  species TEXT,
  tissue TEXT,
  cell_type TEXT,
  dose NUMERIC,
  dose_unit TEXT,
  time_value NUMERIC,
  time_unit TEXT,
  assay TEXT,
  direction TEXT,
  effect_size NUMERIC,
  p_value NUMERIC,
  adjusted_p NUMERIC,
  confidence NUMERIC,
  uncertainty JSONB NOT NULL DEFAULT '{}'::jsonb,
  version TEXT,
  payload JSONB NOT NULL DEFAULT '{}'::jsonb,
  retrieved_at TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX IF NOT EXISTS evidence_edge_subject_idx ON evidence_edge(subject_id, predicate);
CREATE INDEX IF NOT EXISTS evidence_edge_source_idx ON evidence_edge(source_database, source_id);
CREATE TABLE IF NOT EXISTS model_run (
  run_id TEXT PRIMARY KEY,
  task TEXT NOT NULL,
  module TEXT NOT NULL,
  model_name TEXT,
  model_version TEXT,
  status TEXT NOT NULL,
  input_snapshot JSONB NOT NULL DEFAULT '{}'::jsonb,
  parameters JSONB NOT NULL DEFAULT '{}'::jsonb,
  metrics JSONB NOT NULL DEFAULT '{}'::jsonb,
  assumptions JSONB NOT NULL DEFAULT '{}'::jsonb,
  started_at TIMESTAMPTZ,
  finished_at TIMESTAMPTZ,
  code_version TEXT
);
CREATE TABLE IF NOT EXISTS artifact (
  artifact_id TEXT PRIMARY KEY,
  run_id TEXT REFERENCES model_run(run_id),
  artifact_type TEXT NOT NULL,
  uri TEXT NOT NULL,
  sha256 TEXT,
  source_license TEXT,
  metadata JSONB NOT NULL DEFAULT '{}'::jsonb,
  created_at TIMESTAMPTZ NOT NULL DEFAULT now()
);
