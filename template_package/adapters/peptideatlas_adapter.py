"""
PeptideAtlas Adapter for BioCypher.

Loads PeptideAtlas human proteomics data and generates:
- ProteomicsPeptide edges (protein → detected peptide with mass spec evidence)

PeptideAtlas compiles and standardizes peptide identifications from
mass spectrometry proteomics experiments.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class PeptideAtlasAdapter:
    def __init__(self, data_dir="template_package/data/peptideatlas"):
        self.data_dir = Path(data_dir)
        self.proteins = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PeptideAtlas human summary data."""
        path = self.data_dir / 'peptideatlas_human_summary.tsv'
        if not path.exists():
            logger.warning("PeptideAtlas: summary file not found")
            return

        logger.info("PeptideAtlas: Loading proteomics data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = reader.fieldnames or []

            for row in reader:
                # Get whatever fields are available
                protein_id = ''
                for key in ['biosequence_name', 'protein_id', 'accession', 'Entry']:
                    if key in row and row[key].strip():
                        protein_id = row[key].strip()
                        break

                if not protein_id:
                    # Try first column
                    first_key = fieldnames[0] if fieldnames else ''
                    if first_key and row.get(first_key, '').strip():
                        protein_id = row[first_key].strip()

                if not protein_id:
                    continue

                # Collect available data
                data = {}
                for key in ['n_observations', 'n_distinct_peptides', 'protein_group_probability',
                           'presence_level', 'subsumed_by', 'estimated_ng_per_ml',
                           'norm_PSMs_per_100K', 'chromosome']:
                    if key in row:
                        data[key] = row[key].strip()

                self.proteins.append({
                    'id': protein_id,
                    'data': data,
                })
                count += 1

        logger.info(f"PeptideAtlas: Loaded {count} protein entries")

    def get_nodes(self):
        """No new nodes - links to existing proteins."""
        logger.info("PeptideAtlas: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate ProteomicsDetection edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PeptideAtlas: Generating edges...")
        count = 0

        for prot in self.proteins:
            n_obs = prot['data'].get('n_observations', '0')
            n_pep = prot['data'].get('n_distinct_peptides', '0')
            level = prot['data'].get('presence_level', '')

            try:
                n_obs_int = int(n_obs)
            except ValueError:
                n_obs_int = 0

            try:
                n_pep_int = int(n_pep)
            except ValueError:
                n_pep_int = 0

            props = {
                'n_observations': n_obs_int,
                'n_distinct_peptides': n_pep_int,
                'presence_level': level,
                'source': 'PeptideAtlas',
            }

            yield (
                None,
                prot['id'],
                "PEPTIDEATLAS_HUMAN",
                "ProteomicsDetection",
                props
            )
            count += 1

        logger.info(f"PeptideAtlas: Generated {count} ProteomicsDetection edges")
