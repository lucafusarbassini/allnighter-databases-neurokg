"""
OmniPath Adapter for BioCypher.

Loads OmniPath signaling network interactions and generates:
- SignalingInteraction edges (protein → protein signaling)

OmniPath is a comprehensive collection of literature-curated
signaling pathway resources, integrating data from >100 databases.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class OmniPathAdapter:
    def __init__(self, data_dir="template_package/data/omnipath"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load OmniPath signaling interactions."""
        path = self.data_dir / 'omnipath_interactions.tsv'
        if not path.exists():
            logger.warning("OmniPath: interactions file not found")
            return

        logger.info("OmniPath: Loading signaling interactions...")
        count = 0
        seen = set()

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                source = (row.get('source') or '').strip()
                target = (row.get('target') or '').strip()
                source_gene = (row.get('source_genesymbol') or '').strip()
                target_gene = (row.get('target_genesymbol') or '').strip()
                is_directed = (row.get('is_directed') or '').strip()
                is_stimulation = (row.get('is_stimulation') or '').strip()
                is_inhibition = (row.get('is_inhibition') or '').strip()
                sources = (row.get('sources') or '').strip()

                if not source or not target:
                    continue

                key = (source, target)
                if key in seen:
                    continue
                seen.add(key)

                self.interactions.append({
                    'source': source,
                    'target': target,
                    'source_gene': source_gene,
                    'target_gene': target_gene,
                    'is_directed': is_directed,
                    'is_stimulation': is_stimulation,
                    'is_inhibition': is_inhibition,
                    'sources': sources,
                })
                count += 1

        logger.info(f"OmniPath: Loaded {count} signaling interactions")

    def get_nodes(self):
        """No new nodes."""
        logger.info("OmniPath: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate SignalingInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("OmniPath: Generating edges...")
        count = 0

        for ix in self.interactions:
            props = {
                'source_gene': self._sanitize(ix['source_gene']),
                'target_gene': self._sanitize(ix['target_gene']),
                'is_directed': ix['is_directed'],
                'is_stimulation': ix['is_stimulation'],
                'is_inhibition': ix['is_inhibition'],
                'databases': self._sanitize(ix['sources'][:200]),
                'source': 'OmniPath',
            }

            yield (
                None,
                ix['source'],
                ix['target'],
                "SignalingInteraction",
                props
            )
            count += 1

        logger.info(f"OmniPath: Generated {count} SignalingInteraction edges")
