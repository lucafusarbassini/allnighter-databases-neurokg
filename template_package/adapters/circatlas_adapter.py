"""
circAtlas Adapter for BioCypher.

Loads circAtlas v3.0 human circular RNA data and generates:
- CircRNAAtlas nodes (circular RNA annotations with genomic coordinates)

circAtlas is a comprehensive database of circular RNAs across
vertebrate species with tissue expression profiles.
"""

from pathlib import Path
from biocypher._logger import logger


class CircAtlasAdapter:
    def __init__(self, data_dir="template_package/data/circatlas"):
        self.data_dir = Path(data_dir)
        self.circrnas = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load circAtlas human circRNA data."""
        path = self.data_dir / 'human_bed_v3.0.txt'
        if not path.exists():
            logger.warning("circAtlas: BED data not found")
            return

        logger.info("circAtlas: Loading human circRNA data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 5:
                    continue

                chrom = parts[0].strip()
                start = parts[1].strip()
                end = parts[2].strip()
                strand = parts[3].strip()
                circ_id = parts[4].strip()

                if not chrom.startswith('chr'):
                    continue

                try:
                    start_int = int(start)
                    end_int = int(end)
                except ValueError:
                    continue

                self.circrnas.append({
                    'circ_id': circ_id,
                    'chromosome': chrom,
                    'start': start_int,
                    'end': end_int,
                    'strand': strand,
                })
                count += 1

                if count >= 500000:
                    break

        logger.info(f"circAtlas: Loaded {count} human circRNAs")

    def get_nodes(self):
        """
        Generate CircRNAAtlas nodes.
        Yields: (id, label, properties)
        """
        logger.info("circAtlas: Generating nodes...")
        count = 0

        for circ in self.circrnas:
            props = {
                'chromosome': circ['chromosome'],
                'start': circ['start'],
                'end': circ['end'],
                'strand': circ['strand'],
                'source': 'circAtlas_v3',
            }

            yield (f"circAtlas:{circ['circ_id']}", "CircRNAAtlas", props)
            count += 1

        logger.info(f"circAtlas: Generated {count} CircRNAAtlas nodes")

    def get_edges(self):
        """No edges."""
        logger.info("circAtlas: No edges to generate")
        return iter([])
