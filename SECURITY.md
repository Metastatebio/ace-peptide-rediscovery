# Security notes

The lightweight release verifier uses only the Python standard library. The optional analysis
environment fetches AHTPDB from its documented upstream URL and may load model weights only from
their official upstream provider.

Dependency pins are reviewed at release time. `requests` is pinned at a patched version. The
historical ProteinLM requirements deliberately retain the exact `transformers` version recorded for
the frozen analysis; GitHub may flag that reproducibility-only file. Do not use it for new or
network-connected analysis. For new ProteinLM work, use `requirements-proteinlm-maintained.txt`
and revalidate outputs. GitHub currently reports an upstream Biopython advisory without a patched
version. This package does not use `Bio.Entrez` or parse untrusted XML; Biopython is used only for
local peptide properties and optional local PDB handling. Do not process untrusted structure or XML
files in the optional analysis environment.

Please report a repository security issue privately to can@metastate.bio.
