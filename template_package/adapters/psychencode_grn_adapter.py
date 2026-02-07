"""
PsychENCODE GRN Adapter for BioCypher.

Loads PsychENCODE gene regulatory network (ElasticNet filtered) and generates:
- BrainGRN edges (TF → target gene regulatory interactions in brain)

PsychENCODE provides brain-specific gene regulatory networks derived from
Hi-C, ATAC-seq, and RNA-seq data across brain regions.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class PsychENCODEGRNAdapter:
    def __init__(self, data_dir="template_package/data/psychencode"):
        self.data_dir = Path(data_dir)
        self.edges = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PsychENCODE GRN data (filtered ElasticNet)."""
        path = self.data_dir / 'INT-11_ElasticNet_Filtered_Cutoff_0.1_GRN_1.csv'
        if not path.exists():
            logger.warning("PsychENCODE GRN: ElasticNet GRN file not found")
            return

        logger.info("PsychENCODE GRN: Loading gene regulatory network...")
        count = 0
        seen = set()

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                tf = (row.get('Transcription_Factor') or '').strip()
                target = (row.get('Target_Gene') or '').strip()
                enhancer = (row.get('Enhancer_Region') or '').strip()
                weight = (row.get('Edge_Weight') or '').strip()

                if not tf or not target:
                    continue

                # Deduplicate by TF-target pair (keep first/strongest)
                key = (tf, target)
                if key in seen:
                    continue
                seen.add(key)

                try:
                    weight_float = float(weight)
                except ValueError:
                    weight_float = 0.0

                self.edges.append({
                    'tf': tf,
                    'target': target,
                    'enhancer': enhancer,
                    'weight': weight_float,
                })
                count += 1

                if count >= 1000000:
                    break

        logger.info(f"PsychENCODE GRN: Loaded {count} TF-target regulatory edges")

    def get_nodes(self):
        """No new nodes."""
        logger.info("PsychENCODE GRN: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate BrainGRN edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PsychENCODE GRN: Generating edges...")
        count = 0

        for edge in self.edges:
            props = {
                'enhancer_region': self._sanitize(edge['enhancer']),
                'edge_weight': edge['weight'],
                'source': 'PsychENCODE_GRN',
            }

            yield (
                None,
                edge['tf'],
                edge['target'],
                "BrainGRN",
                props
            )
            count += 1

        logger.info(f"PsychENCODE GRN: Generated {count} BrainGRN edges")
