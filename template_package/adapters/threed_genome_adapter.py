"""
3D Genome Browser Adapter for BioCypher.

Loads TAD (Topologically Associating Domain) data from ENCODE
chromatin conformation experiments and generates:
- TAD nodes (genomic domains that represent 3D chromatin organization)

TADs are self-interacting genomic regions that play key roles in
gene regulation and genome organization.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class ThreeDGenomeAdapter:
    def __init__(self, data_dir="template_package/data/3dgenome"):
        self.data_dir = Path(data_dir)
        self.tads = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load TAD BED files from 3D Genome Browser data."""
        logger.info("3DGenome: Loading TAD data...")
        total_count = 0

        for bed_file in sorted(self.data_dir.glob('*.bed.gz')):
            experiment_id = bed_file.stem.replace('.bed', '')
            count = 0

            try:
                with gzip.open(bed_file, 'rt', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#') or line.startswith('track'):
                            continue

                        parts = line.split('\t')
                        if len(parts) < 4:
                            continue

                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                        score = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0

                        tad_id = f"3DG:{experiment_id}:{chrom}:{start}-{end}"

                        self.tads.append({
                            'id': tad_id,
                            'chromosome': chrom,
                            'start': start,
                            'end': end,
                            'name': name,
                            'score': score,
                            'experiment_id': experiment_id,
                        })
                        count += 1
            except Exception as e:
                logger.warning(f"3DGenome: Error reading {bed_file.name}: {e}")
                continue

            logger.info(f"3DGenome: Loaded {count} TADs from {bed_file.name}")
            total_count += count

        logger.info(f"3DGenome: Loaded {total_count} TADs total")

    def get_nodes(self):
        """
        Generate TAD nodes.
        Yields: (id, label, properties)
        """
        logger.info("3DGenome: Generating TAD nodes...")
        count = 0

        for tad in self.tads:
            props = {
                'chromosome': tad['chromosome'],
                'start': tad['start'],
                'end': tad['end'],
                'name': self._sanitize(tad['name']),
                'score': tad['score'],
                'experiment_id': tad['experiment_id'],
                'genome_assembly': 'GRCh38',
                'source': '3DGenomeBrowser',
            }

            yield (tad['id'], "TopologicalDomain", props)
            count += 1

        logger.info(f"3DGenome: Generated {count} TAD nodes")

    def get_edges(self):
        """No edges for TADs currently."""
        logger.info("3DGenome: No edges to generate")
        return iter([])
