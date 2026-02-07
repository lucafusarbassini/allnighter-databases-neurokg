"""
HistoneDB Adapter for BioCypher.

Loads HistoneDB histone variant data and generates:
- HistoneVariant nodes (histone proteins with variant classification)

HistoneDB catalogs histone protein sequences classified by type
(H1, H2A, H2B, H3, H4) and variant.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class HistoneDBAdapter:
    def __init__(self, data_dir="template_package/data/histonedb"):
        self.data_dir = Path(data_dir)
        self.histones = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load HistoneDB human sequences."""
        path = self.data_dir / 'human_all_sequences_complete.json'
        if not path.exists():
            logger.warning("HistoneDB: data file not found")
            return

        logger.info("HistoneDB: Loading histone data...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        rows = data.get('rows', [])
        count = 0
        seen = set()

        for row in rows:
            seq_id = row.get('id', '').strip()
            variant = row.get('variant', '').strip()
            histone_type = row.get('type', '').strip()
            gene = row.get('gene', '')
            score = row.get('score', 0)

            if not seq_id or seq_id in seen:
                continue
            seen.add(seq_id)

            self.histones.append({
                'id': seq_id,
                'variant': variant,
                'histone_type': histone_type,
                'gene': gene or '',
                'score': score or 0,
            })
            count += 1

        logger.info(f"HistoneDB: Loaded {count} histone entries")

    def get_nodes(self):
        """
        Generate HistoneVariant nodes.
        Yields: (id, label, properties)
        """
        logger.info("HistoneDB: Generating nodes...")
        count = 0

        for h in self.histones:
            props = {
                'variant': h['variant'],
                'histone_type': h['histone_type'],
                'gene': h['gene'],
                'score': h['score'],
                'source': 'HistoneDB',
            }

            yield (f"HISTONE:{h['id']}", "HistoneVariant", props)
            count += 1

        logger.info(f"HistoneDB: Generated {count} HistoneVariant nodes")

    def get_edges(self):
        """No edges."""
        logger.info("HistoneDB: No edges to generate")
        return iter([])
