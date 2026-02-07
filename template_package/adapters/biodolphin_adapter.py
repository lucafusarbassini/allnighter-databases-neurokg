"""
Lipid-protein binding Adapter for BioCypher.
Stub adapter for lipid-protein binding data. Ready to load when data becomes available.
"""

import json
import gzip
from pathlib import Path
from biocypher._logger import logger


class BioDolphinAdapter:
    def __init__(self, data_dir="template_package/data/biodolphin"):
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
            logger.warning("BioDolphinAdapter: data directory not found")
            return
        candidates = (list(self.data_dir.glob("*.tsv")) + list(self.data_dir.glob("*.csv"))
                      + list(self.data_dir.glob("*.txt")) + list(self.data_dir.glob("*.json"))
                      + list(self.data_dir.glob("*.tsv.gz")))
        for fpath in candidates:
            try:
                if fpath.suffix == '.json':
                    with open(fpath, 'r') as f:
                        data = json.load(f)
                    if isinstance(data, list):
                        self.entries.extend(data)
                elif str(fpath).endswith('.gz'):
                    with gzip.open(fpath, 'rt', errors='replace') as f:
                        first = f.readline()
                        if first.startswith('<'): continue
                        f.seek(0)
                        self._parse_tsv(f)
                else:
                    with open(fpath, 'r', errors='replace') as f:
                        first = f.readline()
                        if first.startswith('<'): continue
                        f.seek(0)
                        self._parse_tsv(f)
            except Exception as e:
                logger.warning(f"BioDolphinAdapter: Error reading {fpath}: {e}")
        logger.info(f"BioDolphinAdapter: Loaded {len(self.entries)} entries")

    def _parse_tsv(self, fh):
        header = None
        for line in fh:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            if len(parts) >= 2:
                self.entries.append(dict(zip(header, parts)))

    def get_nodes(self):
        logger.info("BioDolphinAdapter: No dedicated nodes")
        return
        yield

    def get_edges(self):
        logger.info("BioDolphinAdapter: No data loaded yet")
        return
        yield
