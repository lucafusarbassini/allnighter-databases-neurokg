"""
GtRNAdb Adapter for BioCypher.

Loads GtRNAdb human tRNA gene data and generates:
- tRNAGene nodes (tRNA genes with anticodon and genomic location)

GtRNAdb catalogs transfer RNA genes identified by tRNAscan-SE
across genomes, providing anticodon, amino acid, and structure info.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class GtRNAdbAdapter:
    def __init__(self, data_dir="template_package/data/gtrnadb"):
        self.data_dir = Path(data_dir)
        self.trnas = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load GtRNAdb tRNA gene data."""
        path = self.data_dir / 'tRNAs_ucsc.txt.gz'
        if not path.exists():
            logger.warning("GtRNAdb: tRNA data file not found")
            return

        logger.info("GtRNAdb: Loading tRNA gene data...")
        count = 0

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue

                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                name = parts[4]
                score = parts[0]
                strand = parts[6]
                amino_acid = parts[7]
                anticodon = parts[8]
                intron_info = parts[9]
                trna_score = parts[10] if len(parts) > 10 else '0'

                try:
                    start_int = int(start)
                    end_int = int(end)
                    score_val = float(trna_score)
                except ValueError:
                    continue

                self.trnas.append({
                    'name': name,
                    'chromosome': chrom,
                    'start': start_int,
                    'end': end_int,
                    'strand': strand,
                    'amino_acid': amino_acid,
                    'anticodon': anticodon,
                    'intron': intron_info,
                    'score': score_val,
                })
                count += 1

        logger.info(f"GtRNAdb: Loaded {count} tRNA genes")

    def get_nodes(self):
        """
        Generate tRNAGene nodes.
        Yields: (id, label, properties)
        """
        logger.info("GtRNAdb: Generating nodes...")
        count = 0

        for trna in self.trnas:
            props = {
                'chromosome': trna['chromosome'],
                'start': trna['start'],
                'end': trna['end'],
                'strand': trna['strand'],
                'amino_acid': trna['amino_acid'],
                'anticodon': trna['anticodon'],
                'has_intron': 'intron' in trna['intron'].lower() and 'no' not in trna['intron'].lower(),
                'score': trna['score'],
                'source': 'GtRNAdb',
            }

            yield (f"tRNA:{trna['name']}", "tRNAGene", props)
            count += 1

        logger.info(f"GtRNAdb: Generated {count} tRNAGene nodes")

    def get_edges(self):
        """No edges."""
        logger.info("GtRNAdb: No edges to generate")
        return iter([])
