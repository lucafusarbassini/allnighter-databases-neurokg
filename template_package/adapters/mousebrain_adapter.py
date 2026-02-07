"""
MouseBrain.org (Linnarsson Lab Single-Cell Atlas) Adapter for BioCypher.

Loads mouse brain cell type taxonomy data.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MouseBrainAdapter:
    def __init__(self, data_dir="template_package/data/mousebrain"):
        self.data_dir = Path(data_dir)
        self.cell_types = []
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
            logger.warning("MouseBrain: data directory not found")
            return
        # Try JSON files
        for fpath in self.data_dir.glob("*.json"):
            try:
                with open(fpath, 'r') as f:
                    data = json.load(f)
                if isinstance(data, list):
                    self.cell_types.extend(data)
                elif isinstance(data, dict):
                    for k, v in data.items():
                        if isinstance(v, dict):
                            v['id'] = k
                            self.cell_types.append(v)
            except Exception as e:
                logger.warning(f"MouseBrain: Error reading {fpath}: {e}")
        # Try TSV files
        for fpath in list(self.data_dir.glob("*.tsv")) + list(self.data_dir.glob("*.csv")):
            try:
                with open(fpath, 'r', errors='replace') as f:
                    first = f.readline()
                    if first.startswith('<'):
                        continue
                    f.seek(0)
                    self._parse_tsv(f)
            except Exception as e:
                logger.warning(f"MouseBrain: Error reading {fpath}: {e}")
        logger.info(f"MouseBrain: Loaded {len(self.cell_types)} cell types")

    def _parse_tsv(self, fh):
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            if len(parts) < 2:
                continue
            row = dict(zip(header, parts))
            self.cell_types.append(row)

    def get_nodes(self):
        logger.info("MouseBrain: Generating CellType nodes...")
        count = 0
        for ct in self.cell_types:
            ct_id = ct.get('id', ct.get('ClusterName', ct.get('cluster_id', '')))
            if not ct_id:
                continue
            name = ct.get('name', ct.get('Description', ct.get('cell_type', '')))
            props = {
                'name': self._sanitize(name),
                'definition': self._sanitize(ct.get('description', ct.get('region', ''))),
                'synonyms': '',
                'xrefs': '',
                'source': 'MouseBrain.org',
            }
            yield (f"mousebrain:{ct_id}", "CellType", props)
            count += 1
        logger.info(f"MouseBrain: Generated {count} CellType nodes")

    def get_edges(self):
        logger.info("MouseBrain: No edges")
        return
        yield
