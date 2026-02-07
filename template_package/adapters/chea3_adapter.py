"""
ChEA3 (ChIP-X Enrichment Analysis 3) Adapter for BioCypher.

Loads ChEA3 transcription factor target gene sets and generates:
- TFTargetInteraction edges (TF → target gene from ChIP experiments)

ChEA3 integrates TF-target gene sets from ChIP-seq/ChIP-chip experiments,
co-expression data, and ENCODE/Consensus libraries.
"""

from pathlib import Path
from biocypher._logger import logger


class ChEA3Adapter:
    def __init__(self, data_dir="template_package/data/chea3"):
        self.data_dir = Path(data_dir)
        self.tf_targets = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ChEA3 GMT files."""
        # Load the primary ChEA library (ChIP-based)
        gmt_files = [
            ('ChEA_2022.gmt', 'ChEA_ChIP'),
            ('ENCODE_ChEA_Consensus_TFs.gmt', 'ENCODE_Consensus'),
        ]

        total = 0
        for filename, library in gmt_files:
            path = self.data_dir / filename
            if not path.exists():
                continue

            logger.info(f"ChEA3: Loading {filename}...")
            count = 0

            with open(path, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue

                    # GMT format: TF_info \t description \t gene1 \t gene2 ...
                    tf_info = parts[0].strip()
                    targets = [g.strip() for g in parts[2:] if g.strip()]

                    # Extract TF name (first word before space or PMID)
                    tf_name = tf_info.split()[0] if tf_info else ''
                    if not tf_name:
                        continue

                    for target in targets:
                        self.tf_targets.append({
                            'tf': tf_name,
                            'target': target,
                            'library': library,
                            'experiment': tf_info,
                        })
                        count += 1

            total += count
            logger.info(f"ChEA3: Loaded {count} TF-target pairs from {filename}")

        logger.info(f"ChEA3: Total {total} TF-target relationships")

    def get_nodes(self):
        """No new nodes - uses existing Gene nodes."""
        logger.info("ChEA3: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate TFTargetInteraction edges (deduplicated by TF+target+library).
        Yields: (id, source, target, label, properties)
        """
        logger.info("ChEA3: Generating edges...")
        seen = set()
        count = 0

        for record in self.tf_targets:
            key = (record['tf'], record['target'], record['library'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'library': record['library'],
                'experiment': self._sanitize(record['experiment'][:200]),
                'source': 'ChEA3',
            }

            yield (
                None,
                record['tf'],
                record['target'],
                "TFTargetInteraction",
                props
            )
            count += 1

        logger.info(f"ChEA3: Generated {count} TFTargetInteraction edges")
