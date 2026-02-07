"""
Enhancer RNA Adapter for BioCypher.

Loads FANTOM5 enhancer region data from BED format and generates:
- EnhancerRegion nodes (genomic regions identified as enhancers)

Data files:
- fantom5_enhancers.bed.gz  (63K enhancer regions in BED12 format)

Each line represents an enhancer region with genomic coordinates,
expression score, and block structure.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class ERNAbaseAdapter:
    def __init__(self, data_dir="template_package/data/ernabase"):
        self.data_dir = Path(data_dir)
        self.enhancers = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load FANTOM5 enhancer BED data."""
        path = self.data_dir / 'fantom5_enhancers.bed.gz'
        if not path.exists():
            logger.warning("ERNAbase: fantom5_enhancers.bed.gz not found")
            return

        logger.info("ERNAbase: Loading FANTOM5 enhancer regions...")
        count = 0

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.split('\t')
                # BED12 format: chr, start, end, name, score, strand,
                #               thickStart, thickEnd, rgb, blockCount, blockSizes, blockStarts
                if len(parts) < 4:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                score = int(parts[4]) if len(parts) > 4 else 0
                strand = parts[5] if len(parts) > 5 else '.'
                thick_start = int(parts[6]) if len(parts) > 6 else start
                thick_end = int(parts[7]) if len(parts) > 7 else end
                block_count = int(parts[9]) if len(parts) > 9 else 1
                block_sizes = parts[10] if len(parts) > 10 else ''
                block_starts = parts[11] if len(parts) > 11 else ''

                self.enhancers.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': strand,
                    'thick_start': thick_start,
                    'thick_end': thick_end,
                    'block_count': block_count,
                    'block_sizes': block_sizes,
                    'block_starts': block_starts,
                    'length': end - start,
                })
                count += 1

        logger.info(f"ERNAbase: Loaded {count} enhancer regions")

    def get_nodes(self):
        """
        Generate EnhancerRegion nodes from FANTOM5 BED data.

        Yields: (id, label, properties)
          - id: enhancer region name (e.g., chr10:100006233-100006603)
        """
        logger.info("ERNAbase: Generating EnhancerRegion nodes...")
        count = 0

        for enh in self.enhancers:
            node_id = enh['name']
            props = {
                'chr': enh['chrom'],
                'start': enh['start'],
                'end': enh['end'],
                'length': enh['length'],
                'score': enh['score'],
                'strand': enh['strand'],
                'thick_start': enh['thick_start'],
                'thick_end': enh['thick_end'],
                'block_count': enh['block_count'],
                'block_sizes': enh['block_sizes'],
                'block_starts': enh['block_starts'],
                'source': 'FANTOM5',
            }

            yield (node_id, "EnhancerRegion", props)
            count += 1

        logger.info(f"ERNAbase: Generated {count} EnhancerRegion nodes")

    def get_edges(self):
        """No edges - enhancer regions are standalone nodes."""
        logger.info("ERNAbase: No edges")
        return iter([])
