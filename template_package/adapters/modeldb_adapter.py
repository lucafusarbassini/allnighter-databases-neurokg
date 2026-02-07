"""
ModelDB Adapter for BioCypher.

Loads ModelDB computational neuroscience model data and generates:
- ComputationalModel nodes (neuroscience models with cell types and channels)

ModelDB is a curated database of published computational neuroscience
models, cataloging cell types, ion channels, and model concepts.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class ModelDBAdapter:
    def __init__(self, data_dir="template_package/data/modeldb"):
        self.data_dir = Path(data_dir)
        self.models = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ModelDB model catalog."""
        path = self.data_dir / 'model_catalog.json'
        if not path.exists():
            logger.warning("ModelDB: catalog file not found")
            return

        logger.info("ModelDB: Loading computational model data...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        count = 0
        for model in data:
            model_id = model.get('id', '')
            name = model.get('name', '').strip()
            cell_types = model.get('cell_types', [])
            ion_channels = model.get('ion_channels', [])
            concepts = model.get('model_concepts', [])

            if not model_id:
                continue

            self.models.append({
                'id': str(model_id),
                'name': name,
                'cell_types': cell_types,
                'ion_channels': ion_channels,
                'concepts': concepts,
            })
            count += 1

        logger.info(f"ModelDB: Loaded {count} computational models")

    def get_nodes(self):
        """
        Generate ComputationalModel nodes.
        Yields: (id, label, properties)
        """
        logger.info("ModelDB: Generating nodes...")
        count = 0

        for model in self.models:
            props = {
                'name': self._sanitize(model['name'][:200]),
                'cell_types': '|'.join(self._sanitize(ct) for ct in model['cell_types'][:10]),
                'ion_channels': '|'.join(self._sanitize(ic) for ic in model['ion_channels'][:10]),
                'model_concepts': '|'.join(self._sanitize(c) for c in model['concepts'][:10]),
                'source': 'ModelDB',
            }

            yield (f"MODELDB:{model['id']}", "ComputationalModel", props)
            count += 1

        logger.info(f"ModelDB: Generated {count} ComputationalModel nodes")

    def get_edges(self):
        """No edges."""
        logger.info("ModelDB: No edges to generate")
        return iter([])
