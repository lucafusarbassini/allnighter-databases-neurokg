"""
JASPAR Adapter for BioCypher.

Loads JASPAR CORE vertebrate transcription factor binding profiles and generates:
- TFBindingProfile nodes (transcription factor binding matrices)

JASPAR is the open-access database of curated, non-redundant transcription
factor (TF) binding profiles stored as position frequency matrices (PFMs).
"""

import json
from pathlib import Path
from biocypher._logger import logger


class JASPARAdapter:
    def __init__(self, data_dir="template_package/data/jaspar"):
        self.data_dir = Path(data_dir)
        self.profiles = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load JASPAR TF binding profiles."""
        path = self.data_dir / 'jaspar_core_vertebrates.json'
        if not path.exists():
            logger.warning("JASPAR: profile data not found")
            return

        logger.info("JASPAR: Loading TF binding profiles...")

        all_results = []
        for json_file in sorted(self.data_dir.glob('jaspar_core_vertebrates*.json')):
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            all_results.extend(data.get('results', []))

        results = all_results
        seen = set()

        for entry in results:
            matrix_id = entry.get('matrix_id', '')
            name = entry.get('name', '')
            base_id = entry.get('base_id', '')
            version = entry.get('version', '')
            collection = entry.get('collection', '')

            if not matrix_id or matrix_id in seen:
                continue
            seen.add(matrix_id)

            self.profiles.append({
                'matrix_id': matrix_id,
                'name': name,
                'base_id': base_id,
                'version': version,
                'collection': collection,
            })

        logger.info(f"JASPAR: Loaded {len(self.profiles)} TF binding profiles")

    def get_nodes(self):
        """
        Generate TFBindingProfile nodes.
        Yields: (id, label, properties)
        """
        logger.info("JASPAR: Generating nodes...")
        count = 0

        for profile in self.profiles:
            props = {
                'name': self._sanitize(profile['name']),
                'base_id': profile['base_id'],
                'version': profile['version'],
                'collection': profile['collection'],
                'source': 'JASPAR',
            }

            yield (f"JASPAR:{profile['matrix_id']}", "TFBindingProfile", props)
            count += 1

        logger.info(f"JASPAR: Generated {count} TFBindingProfile nodes")

    def get_edges(self):
        """No edges."""
        logger.info("JASPAR: No edges to generate")
        return iter([])
