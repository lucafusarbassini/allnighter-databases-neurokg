"""
TISSUES Adapter for BioCypher.

Loads JensenLab TISSUES tissue expression data and generates:
- TissueExpression edges (protein → tissue expression associations)

TISSUES integrates evidence on tissue expression from experiments,
text mining, knowledge bases, and predictions.
"""

from pathlib import Path
from biocypher._logger import logger


class TISSUESAdapter:
    def __init__(self, data_dir="template_package/data/tissues"):
        self.data_dir = Path(data_dir)
        self.associations = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load TISSUES expression data."""
        path = self.data_dir / 'human_tissue_knowledge.tsv'
        if not path.exists():
            logger.warning("TISSUES: data file not found")
            return

        logger.info("TISSUES: Loading tissue expression data...")
        count = 0
        seen = set()

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue

                ensp = parts[0].strip()
                gene = parts[1].strip()
                tissue_id = parts[2].strip()
                tissue_name = parts[3].strip()
                evidence_source = parts[4].strip()
                confidence = parts[5].strip()
                score = parts[6].strip()

                if not ensp or not tissue_id:
                    continue

                key = (ensp, tissue_id)
                if key in seen:
                    continue
                seen.add(key)

                try:
                    score_float = float(score)
                except ValueError:
                    score_float = 0.0

                self.associations.append({
                    'ensp': ensp,
                    'gene': gene,
                    'tissue_id': tissue_id,
                    'tissue_name': tissue_name,
                    'evidence_source': evidence_source,
                    'confidence': confidence,
                    'score': score_float,
                })
                count += 1

        logger.info(f"TISSUES: Loaded {count} protein-tissue associations")

    def get_nodes(self):
        """No new nodes."""
        logger.info("TISSUES: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate TissueExpression edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("TISSUES: Generating edges...")
        count = 0

        for assoc in self.associations:
            props = {
                'gene_symbol': self._sanitize(assoc['gene']),
                'tissue_name': self._sanitize(assoc['tissue_name']),
                'evidence_source': self._sanitize(assoc['evidence_source']),
                'confidence': self._sanitize(assoc['confidence']),
                'score': assoc['score'],
                'source': 'TISSUES_JensenLab',
            }

            yield (
                None,
                assoc['ensp'],
                assoc['tissue_id'],
                "TissueExpression",
                props
            )
            count += 1

        logger.info(f"TISSUES: Generated {count} TissueExpression edges")
