"""
APPRIS Adapter for BioCypher.

Loads APPRIS principal isoform annotations and generates:
- PrincipalIsoform edges (gene → transcript principal isoform annotation)

APPRIS annotates splice isoforms for protein-coding genes, identifying
principal functional isoforms based on structure, function, and conservation.
"""

from pathlib import Path
from biocypher._logger import logger


class APPRISAdapter:
    def __init__(self, data_dir="template_package/data/appris"):
        self.data_dir = Path(data_dir)
        self.isoforms = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load APPRIS principal isoform data."""
        path = self.data_dir / 'appris_principal.txt'
        if not path.exists():
            logger.warning("APPRIS: principal isoform data not found")
            return

        logger.info("APPRIS: Loading principal isoform annotations...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                gene_name = parts[0].strip()
                gene_id = parts[1].strip()
                transcript_id = parts[2].strip()
                ccds_id = parts[3].strip()
                appris_annotation = parts[4].strip()
                mane = parts[5].strip() if len(parts) > 5 else ''

                if not gene_id or not transcript_id:
                    continue

                self.isoforms.append({
                    'gene_name': gene_name,
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                    'ccds_id': ccds_id,
                    'annotation': appris_annotation,
                    'mane': mane,
                })
                count += 1

        logger.info(f"APPRIS: Loaded {count} isoform annotations")

    def get_nodes(self):
        """No new nodes."""
        logger.info("APPRIS: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate PrincipalIsoform edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("APPRIS: Generating edges...")
        count = 0

        for iso in self.isoforms:
            props = {
                'gene_name': self._sanitize(iso['gene_name']),
                'ccds_id': iso['ccds_id'],
                'appris_annotation': iso['annotation'],
                'mane_status': iso['mane'],
                'source': 'APPRIS',
            }

            yield (
                None,
                iso['gene_id'],
                iso['transcript_id'],
                "PrincipalIsoform",
                props
            )
            count += 1

        logger.info(f"APPRIS: Generated {count} PrincipalIsoform edges")
