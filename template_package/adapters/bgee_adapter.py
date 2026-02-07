"""
Bgee Adapter for BioCypher.

Loads Bgee gene expression calls for human and generates:
- GeneExpression edges (gene → anatomical entity expression)

Bgee provides gene expression data across multiple animal species,
integrating RNA-seq, Affymetrix, in situ hybridization, and EST data.
"""

import gzip
import csv
from pathlib import Path
from biocypher._logger import logger


class BgeeAdapter:
    def __init__(self, data_dir="template_package/data/bgee"):
        self.data_dir = Path(data_dir)
        self.expressions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Bgee human expression calls."""
        path = self.data_dir / 'bgee_human_expr.tsv.gz'
        if not path.exists():
            logger.warning("Bgee: data file not found")
            return

        logger.info("Bgee: Loading human gene expression calls...")
        count = 0
        seen = set()
        max_entries = 500000  # Cap for memory management

        try:
            with gzip.open(path, 'rt', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene_id = (row.get('Gene ID') or '').strip()
                    gene_name = (row.get('Gene name') or '').strip().strip('"')
                    anat_id = (row.get('Anatomical entity ID') or '').strip()
                    anat_name = (row.get('Anatomical entity name') or '').strip().strip('"')
                    expression = (row.get('Expression') or '').strip()
                    quality = (row.get('Call quality') or '').strip()
                    score = (row.get('Expression score') or '').strip()

                    if not gene_id or not anat_id:
                        continue

                    # Only include present expression calls with gold quality
                    if expression != 'present':
                        continue
                    if quality != 'gold quality':
                        continue

                    key = (gene_id, anat_id)
                    if key in seen:
                        continue
                    seen.add(key)

                    try:
                        score_float = float(score)
                    except ValueError:
                        score_float = 0.0

                    self.expressions.append({
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'anat_id': anat_id,
                        'anat_name': anat_name,
                        'score': score_float,
                    })
                    count += 1

                    if count >= max_entries:
                        break

        except EOFError:
            logger.warning(f"Bgee: Truncated gzip, loaded {count} entries")

        logger.info(f"Bgee: Loaded {count} gold-quality expression calls")

    def get_nodes(self):
        """No new nodes."""
        logger.info("Bgee: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate GeneExpression edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("Bgee: Generating edges...")
        count = 0

        for expr in self.expressions:
            props = {
                'gene_name': self._sanitize(expr['gene_name']),
                'anatomical_entity': self._sanitize(expr['anat_name']),
                'expression_score': expr['score'],
                'source': 'Bgee',
            }

            yield (
                None,
                expr['gene_id'],
                expr['anat_id'],
                "GeneExpression",
                props
            )
            count += 1

        logger.info(f"Bgee: Generated {count} GeneExpression edges")
