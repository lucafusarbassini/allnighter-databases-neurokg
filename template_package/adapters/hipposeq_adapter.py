"""
HippoSeq Adapter for BioCypher.
Stub adapter for HippoSeq hippocampus RNA-seq atlas data.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class HippoSeqAdapter:
    def __init__(self, data_dir="template_package/data/hipposeq"):
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
            logger.warning("HippoSeqAdapter: data directory not found")
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
                logger.warning(f"HippoSeqAdapter: Error reading {fpath}: {e}")
        logger.info(f"HippoSeqAdapter: Loaded {len(self.entries)} entries")

    def get_nodes(self):
        logger.info("HippoSeqAdapter: No data loaded (data not available for download)")
        return
        yield

    def get_edges(self):
        logger.info("HippoSeqAdapter: No edges")
        return
        yield
