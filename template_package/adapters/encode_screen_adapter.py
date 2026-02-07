"""
ENCODE SCREEN (cCRE) Adapter for BioCypher.

Loads candidate cis-Regulatory Elements from ENCODE SCREEN and generates:
- CisRegulatoryElement nodes (genomic cCREs with chromosomal coordinates)

The BED file contains ~1M cCREs. For scalability, we load all of them
as they represent the regulatory landscape of the human genome.
"""

from pathlib import Path
from biocypher._logger import logger


class ENCODESCREENAdapter:
    def __init__(self, data_dir="template_package/data/encode_screen"):
        self.data_dir = Path(data_dir)
        self.ccres = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load GRCh38-cCREs.bed file."""
        bed_path = self.data_dir / 'GRCh38-cCREs.bed'
        if not bed_path.exists():
            logger.warning("ENCODE SCREEN: cCRE BED file not found")
            return

        logger.info("ENCODE SCREEN: Loading cCREs...")
        count = 0

        with open(bed_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.split('\t')
                if len(parts) < 6:
                    continue

                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                ccre_id = parts[3]     # EH38D* ID
                screen_id = parts[4]   # EH38E* ID
                ccre_type = parts[5]   # e.g. "pELS,CTCF-bound"

                self.ccres.append({
                    'ccre_id': ccre_id,
                    'screen_id': screen_id,
                    'chromosome': chrom,
                    'start': int(start),
                    'end': int(end),
                    'ccre_type': ccre_type,
                })
                count += 1

        logger.info(f"ENCODE SCREEN: Loaded {count} cCREs")

    def get_nodes(self):
        """
        Generate CisRegulatoryElement nodes.
        Yields: (id, label, properties)
        """
        logger.info("ENCODE SCREEN: Generating nodes...")
        count = 0

        for ccre in self.ccres:
            props = {
                'screen_id': ccre['screen_id'],
                'chromosome': ccre['chromosome'],
                'start': ccre['start'],
                'end': ccre['end'],
                'ccre_type': ccre['ccre_type'],
                'genome_assembly': 'GRCh38',
                'source': 'ENCODE_SCREEN',
            }

            yield (ccre['ccre_id'], "CisRegulatoryElement", props)
            count += 1

        logger.info(f"ENCODE SCREEN: Generated {count} CisRegulatoryElement nodes")

    def get_edges(self):
        """
        No edges for now - cCRE-gene links require additional data.
        Yields nothing.
        """
        logger.info("ENCODE SCREEN: No edges to generate")
        return iter([])
