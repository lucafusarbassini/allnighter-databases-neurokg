"""
BrainKB Adapter for BioCypher.
Stub adapter for BrainKB neuroscience knowledge base. Requires authentication to access.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class BrainKBAdapter:
    def __init__(self, data_dir="template_package/data/brainkb"):
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
            logger.warning("BrainKBAdapter: data directory not found")
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
                logger.warning(f"BrainKBAdapter: Error reading {fpath}: {e}")
        logger.info(f"BrainKBAdapter: Loaded {len(self.entries)} entries")

    def get_nodes(self):
        logger.info("BrainKBAdapter: No data loaded (authentication required)")
        return
        yield

    def get_edges(self):
        logger.info("BrainKBAdapter: No edges")
        return
        yield
