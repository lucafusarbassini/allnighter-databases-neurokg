"""
GO-CAM (Gene Ontology Causal Activity Models) Adapter for BioCypher.

Loads structured causal models from GO-CAM.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class GOCAMAdapter:
    def __init__(self, data_dir="template_package/data/gocam"):
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
        if not self.data_dir.exists():
            logger.warning("GO-CAM: data directory not found")
            return
        for fpath in self.data_dir.glob("*.json"):
            try:
                with open(fpath, 'r') as f:
                    first = f.readline()
                    if first.startswith('<'):
                        continue
                    f.seek(0)
                    data = json.load(f)
                if isinstance(data, list):
                    self.models.extend(data)
                elif isinstance(data, dict):
                    # Handle SPARQL results format
                    if 'results' in data and 'bindings' in data['results']:
                        for binding in data['results']['bindings']:
                            model_uri = binding.get('model', {}).get('value', '')
                            title = binding.get('title', {}).get('value', '')
                            if model_uri and 'go-graphstore' not in model_uri:
                                model_id = model_uri.rstrip('/').split('/')[-1]
                                self.models.append({
                                    'id': model_id,
                                    'title': title,
                                    'uri': model_uri,
                                })
                    elif 'models' in data:
                        self.models.extend(data['models'])
                    else:
                        self.models.append(data)
            except Exception as e:
                logger.warning(f"GO-CAM: Error reading {fpath}: {e}")
        logger.info(f"GO-CAM: Loaded {len(self.models)} models")

    def get_nodes(self):
        logger.info("GO-CAM: Generating CausalActivityModel nodes...")
        count = 0
        for model in self.models:
            model_id = model.get('id', model.get('gocam', ''))
            if not model_id:
                continue
            title = model.get('title', model.get('label', ''))
            props = {
                'name': self._sanitize(title),
                'state': self._sanitize(model.get('state', '')),
                'source': 'GO-CAM',
            }
            yield (model_id, "CausalActivityModel", props)
            count += 1
        logger.info(f"GO-CAM: Generated {count} CausalActivityModel nodes")

    def get_edges(self):
        logger.info("GO-CAM: No edges (models are standalone)")
        return
        yield
