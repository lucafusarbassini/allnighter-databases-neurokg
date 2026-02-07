"""
EpiMap Adapter for BioCypher.

Loads EpiMap (Epigenome Integration across Multiple Annotation Projects)
chromatin state segmentation data. EpiMap extends the Roadmap Epigenomics
project with imputed chromatin state data across hundreds of biosamples.

Generates:
- ChromatinState nodes (genomic regions annotated with chromatin states)

Chromatin states include: TssA, TssFlnk, EnhA1, EnhA2, EnhG1, EnhG2,
TxWk, Tx, ZNF/Rpts, Het, TssBiv, EnhBiv, ReprPC, ReprPCWk, Quies, etc.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class EpiMapAdapter:
    def __init__(self, data_dir="template_package/data/epimap"):
        self.data_dir = Path(data_dir)
        self.states = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load EpiMap chromatin state segmentation files."""
        logger.info("EpiMap: Loading chromatin state data...")
        total = 0

        # Active chromatin state types (skip quiescent for size)
        active_states = {
            'TssA', 'TssFlnk', 'TssFlnkU', 'TssFlnkD',
            'Tx', 'TxWk', 'EnhG1', 'EnhG2', 'EnhA1', 'EnhA2',
            'EnhWk', 'TssBiv', 'EnhBiv', 'ReprPC', 'ReprPCWk',
            'ZNF/Rpts', 'Het'
        }

        for bed_file in sorted(self.data_dir.glob('*.bed.gz')):
            biosample_id = bed_file.stem.split('_')[0]
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
                        state = parts[3]

                        # Only keep active/functional states (not Quies)
                        if state not in active_states:
                            continue

                        state_id = f"EPI:{biosample_id}:{chrom}:{start}-{end}"

                        self.states.append({
                            'id': state_id,
                            'chromosome': chrom,
                            'start': start,
                            'end': end,
                            'state': state,
                            'biosample_id': biosample_id,
                        })
                        count += 1
            except Exception as e:
                logger.warning(f"EpiMap: Error reading {bed_file.name}: {e}")

            logger.info(f"EpiMap: Loaded {count} active states from {bed_file.name}")
            total += count

        logger.info(f"EpiMap: Loaded {total} chromatin state regions total")

    def get_nodes(self):
        """
        Generate ChromatinState nodes.
        Yields: (id, label, properties)
        """
        logger.info("EpiMap: Generating chromatin state nodes...")
        count = 0

        for state in self.states:
            props = {
                'chromosome': state['chromosome'],
                'start': state['start'],
                'end': state['end'],
                'state': state['state'],
                'biosample_id': state['biosample_id'],
                'genome_assembly': 'GRCh38',
                'source': 'EpiMap',
            }

            yield (state['id'], "ChromatinState", props)
            count += 1

        logger.info(f"EpiMap: Generated {count} chromatin state nodes")

    def get_edges(self):
        """No edges generated currently."""
        logger.info("EpiMap: No edges to generate")
        return iter([])
