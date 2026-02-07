"""
SynGO (Synaptic Gene Ontologies) Adapter for BioCypher.

Loads SynGO curated synaptic gene annotations and generates:
- SynapticTerm nodes (GO terms curated for synaptic biology)
- GeneAnnotatedToSynapticTerm edges (gene → GO term with evidence)
- SynapticTermIsA edges (term hierarchy)
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class SynGOAdapter:
    def __init__(self, data_dir="template_package/data/syngo"):
        self.data_dir = Path(data_dir)
        self.terms = {}        # GO:XXXXXXX -> {name, domain, shortname, parent_id}
        self.annotations = []  # [{uniprot_id, go_id, go_name, evidence_*, ...}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load SynGO TSV files."""
        self._load_ontologies()
        self._load_annotations()

    def _load_ontologies(self):
        """Load syngo_ontologies.tsv for synaptic GO terms."""
        path = self.data_dir / 'syngo_ontologies.tsv'
        if not path.exists():
            logger.warning("SynGO: ontologies file not found")
            return

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                go_id = row.get('id', '').strip()
                if not go_id.startswith('GO:'):
                    continue

                # Extract clean name (remove GO ID suffix like "synapse (GO:0045202)")
                raw_name = row.get('shortname', '') or row.get('name', '')
                name = self._sanitize(raw_name)

                self.terms[go_id] = {
                    'name': name,
                    'domain': row.get('domain', '').strip(),  # CC or BP
                    'parent_id': row.get('parent_id', '').strip(),
                }

        logger.info(f"SynGO: Loaded {len(self.terms)} synaptic GO terms")

    def _load_annotations(self):
        """Load syngo_annotations.tsv for gene-to-GO-term mappings."""
        path = self.data_dir / 'syngo_annotations.tsv'
        if not path.exists():
            logger.warning("SynGO: annotations file not found")
            return

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                uniprot_id = row.get('uniprot_id', '').strip()
                go_id = row.get('go_id', '').strip()

                if not uniprot_id or not go_id:
                    continue

                self.annotations.append({
                    'uniprot_id': uniprot_id,
                    'go_id': go_id,
                    'go_name': self._sanitize(row.get('go_name', '')),
                    'go_domain': row.get('go_domain', '').strip(),
                    'biological_system': self._sanitize(row.get('evidence_biological_system', '')),
                    'protein_targeting': self._sanitize(row.get('evidence_protein_targeting', '')),
                    'experiment_assay': self._sanitize(row.get('evidence_experiment_assay', '')),
                    'pubmed_id': row.get('pubmed_id', '').strip(),
                })

        logger.info(f"SynGO: Loaded {len(self.annotations)} gene annotations")

    def get_nodes(self):
        """
        Generate SynapticTerm nodes.
        Yields: (id, label, properties)
        """
        logger.info("SynGO: Generating nodes...")
        count = 0

        for go_id, data in self.terms.items():
            props = {
                'name': data['name'],
                'domain': data['domain'],
                'source': 'SynGO',
            }
            yield (go_id, "SynapticTerm", props)
            count += 1

        logger.info(f"SynGO: Generated {count} SynapticTerm nodes")

    def get_edges(self):
        """
        Generate edges:
        - GeneAnnotatedToSynapticTerm (gene → GO term)
        - SynapticTermIsA (term hierarchy)
        Yields: (id, source, target, label, properties)
        """
        logger.info("SynGO: Generating edges...")

        # Gene-to-term annotations (deduplicated by gene+term)
        seen = set()
        anno_count = 0
        for anno in self.annotations:
            key = (anno['uniprot_id'], anno['go_id'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'go_domain': anno['go_domain'],
                'biological_system': anno['biological_system'],
                'experiment_assay': anno['experiment_assay'],
                'pubmed_id': anno['pubmed_id'],
                'source': 'SynGO',
            }

            yield (
                None,
                anno['uniprot_id'],
                anno['go_id'],
                "GeneAnnotatedToSynapticTerm",
                props
            )
            anno_count += 1

        logger.info(f"SynGO: Generated {anno_count} GeneAnnotatedToSynapticTerm edges")

        # Term hierarchy
        hier_count = 0
        for go_id, data in self.terms.items():
            parent_id = data.get('parent_id', '')
            if parent_id and parent_id.startswith('GO:'):
                yield (None, go_id, parent_id, "SynapticTermIsA", {})
                hier_count += 1

        logger.info(f"SynGO: Generated {hier_count} SynapticTermIsA edges")
