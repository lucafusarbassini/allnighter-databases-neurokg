"""
ConsensusPathDB (CPDB) Adapter for BioCypher.

Loads ConsensusPathDB human protein-protein interactions and generates:
- ConsensusInteraction edges (protein-protein interactions from multiple sources)

CPDB integrates interactions from 32 databases including Reactome, KEGG,
HPRD, IntAct, BioGRID, and others, providing confidence-scored interactions.
"""

import gzip
import csv
from pathlib import Path
from biocypher._logger import logger


class CPDBAdapter:
    def __init__(self, data_dir="template_package/data/cpdb"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load CPDB human PPI data."""
        ppi_path = self.data_dir / 'ConsensusPathDB_human_PPI.gz'
        if not ppi_path.exists():
            logger.warning("CPDB: PPI file not found")
            return

        logger.info("CPDB: Loading human protein-protein interactions...")
        count = 0
        skipped = 0

        with gzip.open(ppi_path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    skipped += 1
                    continue

                source_dbs = parts[0].strip()
                publications = parts[1].strip()
                uniprot_entries = parts[2].strip()  # Entry names like MK03_HUMAN
                uniprot_ids = parts[3].strip()       # Accessions like P27361
                gene_names = parts[4].strip()
                confidence = parts[8].strip() if len(parts) > 8 else ''

                # Parse participants (comma-separated)
                ids = uniprot_ids.split(',')
                genes = gene_names.split(',')

                # Only keep binary interactions (2 participants)
                if len(ids) != 2:
                    skipped += 1
                    continue

                id_a = ids[0].strip()
                id_b = ids[1].strip()

                if not id_a or not id_b:
                    continue

                self.interactions.append({
                    'id_a': id_a,
                    'id_b': id_b,
                    'source_dbs': source_dbs,
                    'num_publications': len(publications.split(',')) if publications else 0,
                    'confidence': confidence,
                })
                count += 1

        logger.info(f"CPDB: Loaded {count} binary interactions (skipped {skipped} non-binary)")

    def get_nodes(self):
        """
        No new nodes - CPDB uses UniProt IDs linking to existing Gene nodes.
        """
        logger.info("CPDB: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate ConsensusInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("CPDB: Generating edges...")
        seen = set()
        count = 0

        for inter in self.interactions:
            # Deduplicate: A-B same as B-A
            key = tuple(sorted([inter['id_a'], inter['id_b']]))
            if key in seen:
                continue
            seen.add(key)

            props = {
                'source_databases': self._sanitize(inter['source_dbs']),
                'num_publications': inter['num_publications'],
                'confidence': inter['confidence'],
                'source': 'CPDB',
            }

            yield (
                None,
                inter['id_a'],
                inter['id_b'],
                "ConsensusInteraction",
                props
            )
            count += 1

        logger.info(f"CPDB: Generated {count} ConsensusInteraction edges")
