"""
Mouse Phenome Database Adapter for BioCypher.
Stub adapter for Mouse Phenome Database phenotypic data.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MousePhenomeAdapter:
    def __init__(self, data_dir="template_package/data/mouse_phenome"):
        self.data_dir = Path(data_dir)
        self.entries = []
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
            logger.warning("MousePhenomeAdapter: data directory not found")
            return
        for fpath in self.data_dir.glob("*.json"):
            try:
                with open(fpath, 'r') as f:
                    data = json.load(f)
                if isinstance(data, list):
                    self.entries.extend(data)
                elif isinstance(data, dict) and 'results' in data:
                    self.entries.extend(data['results'])
            except Exception as e:
                logger.warning(f"MousePhenomeAdapter: Error reading {fpath}: {e}")
        logger.info(f"MousePhenomeAdapter: Loaded {len(self.entries)} entries")

    def get_nodes(self):
        logger.info("MousePhenomeAdapter: No data loaded (data not available for download)")
        return
        yield

    def get_edges(self):
        logger.info("MousePhenomeAdapter: No edges")
        return
        yield
