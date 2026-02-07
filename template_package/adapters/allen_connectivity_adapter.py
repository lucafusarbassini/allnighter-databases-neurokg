"""
Allen Brain Connectivity Atlas Adapter for BioCypher.

Loads Allen Institute mouse brain connectivity experiment data and generates:
- ConnectivityExperiment nodes (tract-tracing experiments)

The Allen Mouse Brain Connectivity Atlas maps neural connections
through viral tracer injections and whole-brain serial two-photon imaging.
"""

import csv
import json
from pathlib import Path
from biocypher._logger import logger


class AllenConnectivityAdapter:
    def __init__(self, data_dir="template_package/data/allen_connectivity"):
        self.data_dir = Path(data_dir)
        self.experiments = []
        self.projections = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Allen Connectivity data."""
        # Load experiment summary
        csv_path = self.data_dir / 'experiment_summary.csv'
        if csv_path.exists():
            logger.info("Allen Connectivity: Loading experiment summary...")
            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    exp_id = row.get('experiment_id', '').strip()
                    failed = row.get('failed', '').strip()

                    if not exp_id or failed == 'True':
                        continue

                    self.experiments.append({
                        'id': exp_id,
                        'specimen_id': row.get('specimen_id', '').strip(),
                        'section_thickness': row.get('section_thickness', '').strip(),
                    })

            logger.info(f"Allen Connectivity: Loaded {len(self.experiments)} experiments")

        # Load projection summary
        proj_path = self.data_dir / 'projection_summary.csv'
        if proj_path.exists():
            logger.info("Allen Connectivity: Loading projection summary...")
            with open(proj_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    self.projections.append(row)

            logger.info(f"Allen Connectivity: Loaded {len(self.projections)} projection entries")

    def get_nodes(self):
        """
        Generate ConnectivityExperiment nodes.
        Yields: (id, label, properties)
        """
        logger.info("Allen Connectivity: Generating nodes...")
        count = 0

        for exp in self.experiments:
            props = {
                'specimen_id': exp['specimen_id'],
                'section_thickness': exp['section_thickness'],
                'source': 'Allen_Connectivity',
            }

            yield (f"ALLEN_CONN:{exp['id']}", "ConnectivityExperiment", props)
            count += 1

        logger.info(f"Allen Connectivity: Generated {count} ConnectivityExperiment nodes")

    def get_edges(self):
        """No edges for now."""
        logger.info("Allen Connectivity: No edges to generate")
        return iter([])
