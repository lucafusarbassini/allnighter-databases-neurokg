"""
PolyASite Adapter for BioCypher.

Loads PolyASite 2.0 polyadenylation site clusters and generates:
- PolyASite nodes (polyadenylation site clusters with tissue expression)

PolyASite provides a comprehensive atlas of polyadenylation sites
across human tissues, mapped from 3' end sequencing data.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class PolyASiteAdapter:
    def __init__(self, data_dir="template_package/data/polyasite"):
        self.data_dir = Path(data_dir)
        self.sites = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PolyASite cluster data."""
        path = self.data_dir / 'polyasite_human.tsv.gz'
        if not path.exists():
            logger.warning("PolyASite: data file not found")
            return

        logger.info("PolyASite: Loading polyadenylation site clusters...")
        count = 0

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')

                if header is None:
                    header = parts
                    continue

                if len(parts) < 11:
                    continue

                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                name = parts[3]
                score = parts[4]
                strand = parts[5]
                rep_site = parts[6]
                frac_samples = parts[7]
                nr_prots = parts[8]
                annotation = parts[9]
                gene_name = parts[10]

                try:
                    start_int = int(start)
                    end_int = int(end)
                    score_float = float(score)
                except ValueError:
                    continue

                self.sites.append({
                    'name': name,
                    'chromosome': chrom,
                    'start': start_int,
                    'end': end_int,
                    'strand': strand,
                    'score': score_float,
                    'representative_site': rep_site,
                    'fraction_samples': frac_samples,
                    'num_protocols': nr_prots,
                    'annotation': annotation,
                    'gene_name': gene_name,
                })
                count += 1

                if count >= 500000:
                    break

        logger.info(f"PolyASite: Loaded {count} polyadenylation site clusters")

    def get_nodes(self):
        """
        Generate PolyASite nodes.
        Yields: (id, label, properties)
        """
        logger.info("PolyASite: Generating nodes...")
        count = 0

        for site in self.sites:
            props = {
                'chromosome': site['chromosome'],
                'start': site['start'],
                'end': site['end'],
                'strand': site['strand'],
                'score': site['score'],
                'annotation': site['annotation'],
                'gene_name': self._sanitize(site['gene_name']),
                'num_protocols': site['num_protocols'],
                'source': 'PolyASite',
            }

            yield (f"PAS:{site['name']}", "PolyACluster", props)
            count += 1

        logger.info(f"PolyASite: Generated {count} PolyACluster nodes")

    def get_edges(self):
        """No edges."""
        logger.info("PolyASite: No edges to generate")
        return iter([])
