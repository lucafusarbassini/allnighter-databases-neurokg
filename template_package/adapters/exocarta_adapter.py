"""
ExoCarta (Exosome Database) Adapter for BioCypher.

Loads ExoCarta exosome cargo data and generates:
- ExosomeProtein edges (gene → exosome detection with species and methods)

ExoCarta catalogs proteins, mRNAs, and lipids identified in exosomes
from various cell types and biological fluids.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class ExoCartaAdapter:
    def __init__(self, data_dir="template_package/data/exocarta"):
        self.data_dir = Path(data_dir)
        self.proteins = []   # [{gene_symbol, species, methods, content_type}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ExoCarta protein/mRNA details."""
        path = self.data_dir / 'EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt'
        if not path.exists():
            logger.warning("ExoCarta: protein/mRNA file not found")
            return

        logger.info("ExoCarta: Loading exosome cargo data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_symbol = row.get('GENE SYMBOL', '').strip()
                species = row.get('SPECIES', '').strip()
                content_type = row.get('CONTENT TYPE', '').strip()
                methods = row.get('METHODS', '').strip()
                entrez_id = row.get('ENTREZ GENE ID', '').strip()

                if not gene_symbol:
                    continue

                # Focus on human and mouse
                if 'sapiens' not in species and 'musculus' not in species:
                    continue

                self.proteins.append({
                    'gene_symbol': gene_symbol,
                    'entrez_id': entrez_id,
                    'species': species,
                    'content_type': content_type,
                    'methods': self._sanitize(methods),
                })
                count += 1

        logger.info(f"ExoCarta: Loaded {count} exosome cargo entries (human/mouse)")

    def get_nodes(self):
        """
        No new nodes - ExoCarta genes link to existing Gene nodes.
        Yields nothing.
        """
        logger.info("ExoCarta: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate GeneInExosome edges (deduplicated by gene+species).
        Yields: (id, source, target, label, properties)
        """
        logger.info("ExoCarta: Generating edges...")

        # Aggregate: for each gene+species, collect content types and methods
        gene_data = {}
        for p in self.proteins:
            key = (p['gene_symbol'], p['species'])
            if key not in gene_data:
                gene_data[key] = {
                    'entrez_id': p['entrez_id'],
                    'content_types': set(),
                    'methods': set(),
                    'count': 0,
                }
            gene_data[key]['content_types'].add(p['content_type'])
            for m in p['methods'].split('|'):
                if m.strip():
                    gene_data[key]['methods'].add(m.strip())
            gene_data[key]['count'] += 1

        count = 0
        for (gene_symbol, species), data in gene_data.items():
            props = {
                'species': species,
                'content_types': '|'.join(sorted(data['content_types'])),
                'detection_methods': '|'.join(sorted(data['methods'])),
                'num_experiments': data['count'],
                'source': 'ExoCarta',
            }

            yield (
                None,
                gene_symbol,
                f"EXOSOME:{species.replace(' ', '_')}",
                "GeneInExosome",
                props
            )
            count += 1

        logger.info(f"ExoCarta: Generated {count} GeneInExosome edges")
