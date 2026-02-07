"""
TarBase / miRDB Adapter for BioCypher.

Loads miRNA target prediction data from miRDB v6.0 (which was
originally intended to complement DIANA-TarBase data).

miRDB provides computationally predicted miRNA targets with
prediction scores. The data includes human miRNA-target pairs.

Generates:
- miRNA nodes
- miRNA-target prediction edges
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class TarBaseAdapter:
    def __init__(self, data_dir="template_package/data/tarbase"):
        self.data_dir = Path(data_dir)
        self.predictions = []
        self.mirnas = set()
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load miRDB prediction data."""
        pred_path = self.data_dir / 'miRDB_v6.0_prediction_result.txt.gz'
        if not pred_path.exists():
            logger.warning("TarBase/miRDB: Prediction file not found")
            return

        logger.info("TarBase/miRDB: Loading prediction data...")
        count = 0
        human_count = 0

        try:
            with gzip.open(pred_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue

                    mirna_id = parts[0].strip()
                    target_id = parts[1].strip()  # RefSeq ID
                    score = float(parts[2].strip())

                    # Only keep human miRNAs (hsa-miR prefix)
                    if not mirna_id.startswith('hsa-'):
                        continue

                    # Only keep high-confidence predictions (score >= 80)
                    if score < 80:
                        continue

                    self.mirnas.add(mirna_id)
                    self.predictions.append({
                        'mirna_id': mirna_id,
                        'target_id': target_id,
                        'score': score,
                    })
                    human_count += 1
                    count += 1

        except Exception as e:
            logger.warning(f"TarBase/miRDB: Error loading data: {e}")

        logger.info(f"TarBase/miRDB: Loaded {human_count} human high-confidence predictions ({len(self.mirnas)} miRNAs)")

    def get_nodes(self):
        """
        Generate miRNA nodes.
        Yields: (id, label, properties)
        """
        logger.info("TarBase/miRDB: Generating miRNA nodes...")
        count = 0

        for mirna_id in sorted(self.mirnas):
            props = {
                'name': mirna_id,
                'source': 'miRDB',
            }
            yield (f"miRDB:{mirna_id}", "MicroRNA", props)
            count += 1

        logger.info(f"TarBase/miRDB: Generated {count} miRNA nodes")

    def get_edges(self):
        """
        Generate miRNA-target prediction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info(f"TarBase/miRDB: Generating edges from {len(self.predictions)} predictions...")
        count = 0

        for pred in self.predictions:
            props = {
                'prediction_score': pred['score'],
                'source': 'miRDB',
            }

            yield (
                None,
                f"miRDB:{pred['mirna_id']}",
                f"RefSeq:{pred['target_id']}",
                "MiRNATargetPrediction",
                props
            )
            count += 1

        logger.info(f"TarBase/miRDB: Generated {count} prediction edges")
