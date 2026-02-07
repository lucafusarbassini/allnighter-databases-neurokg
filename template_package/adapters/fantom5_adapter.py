"""
FANTOM5 Enhancer Adapter for BioCypher.

Loads FANTOM5 enhancer BED data and generates:
- Enhancer nodes (genomic enhancer regions with coordinates)

FANTOM5 identified active enhancers using bidirectional CAGE
transcription across the human genome.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class FANTOM5Adapter:
    def __init__(self, data_dir="template_package/data/fantom5"):
        self.data_dir = Path(data_dir)
        self.enhancers = []
        self._load_data()

    def _load_data(self):
        """Load FANTOM5 enhancer BED file."""
        bed_path = self.data_dir / 'F5.hg38.enhancers.bed.gz'
        if not bed_path.exists():
            logger.warning("FANTOM5: Enhancer BED file not found")
            return

        logger.info("FANTOM5: Loading enhancers...")
        count = 0

        with gzip.open(bed_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.split('\t')
                if len(parts) < 5:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]   # chr:start-end format
                score = int(parts[4]) if parts[4].isdigit() else 0

                self.enhancers.append({
                    'id': name,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'score': score,
                })
                count += 1

        logger.info(f"FANTOM5: Loaded {count} enhancers")

    def get_nodes(self):
        """
        Generate Enhancer nodes.
        Yields: (id, label, properties)
        """
        logger.info("FANTOM5: Generating nodes...")
        count = 0

        for enh in self.enhancers:
            props = {
                'chromosome': enh['chromosome'],
                'start': enh['start'],
                'end': enh['end'],
                'score': enh['score'],
                'genome_assembly': 'GRCh38',
                'source': 'FANTOM5',
            }

            yield (f"F5:{enh['id']}", "Enhancer", props)
            count += 1

        logger.info(f"FANTOM5: Generated {count} Enhancer nodes")

    def get_edges(self):
        """No edges - enhancer-gene links require additional data."""
        logger.info("FANTOM5: No edges to generate")
        return iter([])
