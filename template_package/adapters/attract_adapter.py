"""
ATtRACT Adapter for BioCypher.

Loads ATtRACT RNA-binding protein motif data and generates:
- RBPMotif nodes (RNA-binding protein sequence motifs)

ATtRACT is a database of RNA-binding protein (RBP) motifs, providing
experimentally determined binding preferences.
"""

from pathlib import Path
from biocypher._logger import logger


class ATtRACTAdapter:
    def __init__(self, data_dir="template_package/data/attract"):
        self.data_dir = Path(data_dir)
        self.motifs = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ATtRACT motif data."""
        path = self.data_dir / 'ATtRACT_db.txt'
        if not path.exists():
            logger.warning("ATtRACT: motif data not found")
            return

        logger.info("ATtRACT: Loading RBP motif data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 8:
                    continue

                gene_name = parts[0].strip()
                gene_id = parts[1].strip()
                organism = parts[3].strip()
                motif = parts[4].strip()
                motif_len = parts[5].strip()
                experiment = parts[6].strip()
                database = parts[7].strip()

                if 'Homo_sapiens' not in organism:
                    continue

                if not gene_name or not motif:
                    continue

                self.motifs.append({
                    'gene_name': gene_name,
                    'gene_id': gene_id,
                    'motif': motif,
                    'motif_length': motif_len,
                    'experiment': experiment,
                    'database': database,
                })
                count += 1

        logger.info(f"ATtRACT: Loaded {count} human RBP motifs")

    def get_nodes(self):
        """
        Generate RBPMotif nodes.
        Yields: (id, label, properties)
        """
        logger.info("ATtRACT: Generating nodes...")
        seen = set()
        count = 0

        for motif in self.motifs:
            key = (motif['gene_name'], motif['motif'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'gene_name': self._sanitize(motif['gene_name']),
                'motif_sequence': motif['motif'],
                'motif_length': motif['motif_length'],
                'experiment': self._sanitize(motif['experiment']),
                'database': motif['database'],
                'source': 'ATtRACT',
            }

            yield (
                f"ATtRACT:{motif['gene_name']}:{motif['motif']}",
                "RBPMotif",
                props
            )
            count += 1

        logger.info(f"ATtRACT: Generated {count} RBPMotif nodes")

    def get_edges(self):
        """No edges."""
        logger.info("ATtRACT: No edges to generate")
        return iter([])
