# Memory - Persistent Learnings

## Architecture
- BioCypher 0.10.1 framework with Neo4j 4.4 backend
- Adapters yield tuples: Nodes (id, label, props), Edges (id, source, target, label, props)
- Schema defined in config/schema_config.yaml using BioLink ontology (sentence case)
- Pipeline: adapters -> BioCypher CSV -> Neo4j import
- "Gene" = unified entity for gene+protein, preferred_id=uniprot
- Human-centric: mouse data projected via orthology mapping

## Data Available (Downloaded) - 62 databases
- **ChEBI**: compounds.tsv.gz, names.tsv.gz, relation.tsv.gz, etc (~62K compounds)
- **Rhea**: rhea-directions.tsv (18,343 master rxns), rhea2ec.tsv, etc
- **LIPIDMAPS**: structures_extended.sdf (49,718 lipids), etc
- **KEGG**: lists/, xrefs/, entries/ (19,571 compounds, 12,384 reactions, 584 pathways)
- **LIANA**: Already has working adapter (liana_adapter.py)
- **3DGenome**: ENCODE TAD BED files (4 experiments)
- **Adhesome**: components.csv, interactions.csv (download error - 404)
- **BRENDA**: enzyme.dat (ExPASy format, ~8000 enzymes with human DR cross-refs)
- **ChapNet**: Cytoscape JSON networks + ChaperoneCorrelation.csv
- **ENCORI**: miRNA-target TSV files (3 miRNA families, ~7K records)
- **EpiMap**: Chromatin state BED.gz files (3 biosamples)
- **EVpedia**: browse_origin.csv, circRNAs_anno.csv, longRNAs_anno.csv
- **LipidOntology**: LION-terms.csv, all-LION-lipid-associations.csv
- **RaftProt**: Raftprot.v2.4.txt (1408 entries, space-delimited)
- **SILVA**: tax_slv_ssu_138.2.txt.gz (rRNA taxonomy)
- **TarBase/miRDB**: miRDB_v6.0_prediction_result.txt.gz (multi-species)
- **UniCarbKB**: GlyGen glycan TSV/JSON files (~584 glycans)
- Plus 46 additional databases with existing adapters

## Key Design Decisions
- ChEBI uses chebi_accession as ID (CHEBI:12345)
- Rhea reactions modeled as NODES (not edges) with 4 directional variants per master
- LIPIDMAPS lipids are subclass of chemical_substance
- KEGG compounds, reactions, pathways each get their own node types
- Cross-references stored as JSON string properties
- Sanitize quotes in string fields for CSV safety
- miRNA nodes use MicroRNA label (shared by ENCORI + TarBase/miRDB)
- RaftProt uses mouse-to-human ortholog mapping
- BRENDA reuses Enzyme/EnzymeCatalyzedBy schema entries (same as ExPASy ENZYME)
- EpiMap filters for active chromatin states only (skips Quies to reduce size)
- TarBase/miRDB filters for human miRNAs with score >= 80

## Learnings
- BioCypher handles deduplication automatically when id=None
- Properties not in schema_config are silently ignored
- str[] type for array properties, separated by pipe in CSV
- Use `is_a` in schema to create custom ontology branches
- Docker pipeline: build -> import -> deploy (scripts/build.sh, scripts/import.sh)
- IMPORTANT: Sanitize quotes in string data (replace " with "")
- Adhesome data files contain "404: Not Found" - adapter handles gracefully
- ENCORI .json files are actually TSV with comment headers
- RaftProt uses space-delimited quoted fields

## Schema Stats
- 102 total schema entries
- 65 adapter entries in create_knowledge_graph.py
- 62 data directories
- All data directories now have corresponding adapters
