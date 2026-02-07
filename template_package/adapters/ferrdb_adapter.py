"""
FerrDb Adapter for BioCypher.

Loads ferroptosis gene data and generates:
- GeneInFerroptosis edges (gene → ferroptosis role)

FerrDb catalogs genes involved in ferroptosis (iron-dependent cell death),
classified as drivers, suppressors, markers, or unclassified regulators.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class FerrDbAdapter:
    def __init__(self, data_dir="template_package/data/ferrdb"):
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
        """Load FerrDb gene data."""
        path = self.data_dir / 'ferrdb_all_genes.tsv'
        if not path.exists():
            logger.warning("FerrDb: gene file not found")
            return

        logger.info("FerrDb: Loading ferroptosis gene data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                symbol = row.get('symbol', '').strip()
                name = row.get('name', '').strip()
                gene_type = row.get('genetype', '').strip()
                experiments = row.get('experiments', '0').strip()

                if not symbol:
                    continue

                try:
                    n_exp = int(experiments)
                except ValueError:
                    n_exp = 0

                self.genes.append({
                    'symbol': symbol,
                    'name': name,
                    'gene_type': gene_type,
                    'n_experiments': n_exp,
                })
                count += 1

        logger.info(f"FerrDb: Loaded {count} ferroptosis genes")

    def get_nodes(self):
        """No new nodes."""
        logger.info("FerrDb: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate GeneInFerroptosis edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("FerrDb: Generating edges...")
        count = 0

        for gene in self.genes:
            props = {
                'gene_name': gene['name'],
                'role': gene['gene_type'],
                'n_experiments': gene['n_experiments'],
                'source': 'FerrDb',
            }

            yield (
                None,
                gene['symbol'],
                f"FERROPTOSIS:{gene['gene_type']}",
                "GeneInFerroptosis",
                props
            )
            count += 1

        logger.info(f"FerrDb: Generated {count} GeneInFerroptosis edges")
