"""
RNAgranuleDB Adapter for BioCypher.

Loads RNA granule protein data and generates:
- GeneInRNAGranule edges (gene → RNA granule association with score)

RNAgranuleDB catalogs proteins associated with RNA granules
(stress granules, P-bodies, and other RNP granules).
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class RNAgranuleDBAdapter:
    def __init__(self, data_dir="template_package/data/rnagranuledb"):
        self.data_dir = Path(data_dir)
        self.genes = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load RNA granule gene data."""
        path = self.data_dir / 'rna_granule_genes_scored.tsv'
        if not path.exists():
            logger.warning("RNAgranuleDB: scored gene file not found")
            return

        logger.info("RNAgranuleDB: Loading RNA granule gene data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('GeneName', '').strip()
                score = row.get('Score', '0').strip()
                tier = row.get('Tier', '').strip()

                if not gene_name:
                    continue

                try:
                    score_val = float(score)
                except ValueError:
                    score_val = 0.0

                self.genes.append({
                    'gene_name': gene_name,
                    'score': score_val,
                    'tier': tier,
                })
                count += 1

        logger.info(f"RNAgranuleDB: Loaded {count} RNA granule genes")

    def get_nodes(self):
        """No new nodes."""
        logger.info("RNAgranuleDB: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate GeneInRNAGranule edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("RNAgranuleDB: Generating edges...")
        count = 0

        for gene in self.genes:
            props = {
                'score': gene['score'],
                'tier': gene['tier'],
                'source': 'RNAgranuleDB',
            }

            yield (
                None,
                gene['gene_name'],
                "RNA_GRANULE",
                "GeneInRNAGranule",
                props
            )
            count += 1

        logger.info(f"RNAgranuleDB: Generated {count} GeneInRNAGranule edges")
