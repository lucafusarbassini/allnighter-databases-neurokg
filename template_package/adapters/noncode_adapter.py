"""
NONCODE Adapter for BioCypher.

Loads NONCODE v6 human lncRNA annotations and generates:
- LncRNA nodes (long non-coding RNA annotations with genomic coordinates)

NONCODE is an integrated knowledge database of non-coding RNAs,
providing comprehensive lncRNA annotations.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class NONCODEAdapter:
    def __init__(self, data_dir="template_package/data/noncode"):
        self.data_dir = Path(data_dir)
        self.lncrnas = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load NONCODE lncRNA BED data."""
        path = self.data_dir / 'noncode_lncgene.bed.gz'
        if not path.exists():
            logger.warning("NONCODE: BED file not found")
            return

        logger.info("NONCODE: Loading lncRNA annotations...")
        count = 0

        try:
            with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue

                    chrom = parts[0]
                    start = parts[1]
                    end = parts[2]
                    name = parts[3]
                    score = parts[4]
                    strand = parts[5]

                    if not chrom.startswith('chr'):
                        continue

                    try:
                        start_int = int(start)
                        end_int = int(end)
                    except ValueError:
                        continue

                    self.lncrnas.append({
                        'name': name,
                        'chromosome': chrom,
                        'start': start_int,
                        'end': end_int,
                        'strand': strand,
                    })
                    count += 1
        except EOFError:
            logger.warning("NONCODE: Truncated gzip file, loaded partial data")

        logger.info(f"NONCODE: Loaded {count} lncRNA annotations")

    def get_nodes(self):
        """
        Generate LncRNA nodes.
        Yields: (id, label, properties)
        """
        logger.info("NONCODE: Generating nodes...")
        count = 0

        for lncrna in self.lncrnas:
            props = {
                'chromosome': lncrna['chromosome'],
                'start': lncrna['start'],
                'end': lncrna['end'],
                'strand': lncrna['strand'],
                'source': 'NONCODEv6',
            }

            yield (f"NONCODE:{lncrna['name']}", "LncRNA", props)
            count += 1

        logger.info(f"NONCODE: Generated {count} LncRNA nodes")

    def get_edges(self):
        """No edges."""
        logger.info("NONCODE: No edges to generate")
        return iter([])
