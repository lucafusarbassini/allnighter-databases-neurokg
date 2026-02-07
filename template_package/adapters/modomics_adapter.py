"""
MODOMICS Adapter for BioCypher.

Loads MODOMICS RNA modification data and generates:
- RNAModification nodes (known RNA nucleoside modifications)

MODOMICS is a database of RNA modification pathways, providing
information on modified nucleosides found in RNA.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MODOMICSAdapter:
    def __init__(self, data_dir="template_package/data/modomics"):
        self.data_dir = Path(data_dir)
        self.modifications = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load MODOMICS modification data."""
        path = self.data_dir / 'modifications.json'
        if not path.exists():
            logger.warning("MODOMICS: modification data not found")
            return

        logger.info("MODOMICS: Loading RNA modifications...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Data is a dict keyed by ID
        for mod_id, mod in data.items():
            if isinstance(mod, dict):
                self.modifications.append({
                    'id': mod.get('id', mod_id),
                    'name': mod.get('name', ''),
                    'short_name': mod.get('short_name', ''),
                    'formula': mod.get('formula', ''),
                    'mass_avg': mod.get('mass_avg', 0),
                    'mass_monoiso': mod.get('mass_monoiso', 0),
                    'reference_moiety': ';'.join(mod.get('reference_moiety', [])) if isinstance(mod.get('reference_moiety'), list) else str(mod.get('reference_moiety', '')),
                })

        logger.info(f"MODOMICS: Loaded {len(self.modifications)} RNA modifications")

    def get_nodes(self):
        """
        Generate RNAModification nodes.
        Yields: (id, label, properties)
        """
        logger.info("MODOMICS: Generating nodes...")
        count = 0

        for mod in self.modifications:
            props = {
                'name': self._sanitize(mod['name']),
                'short_name': mod['short_name'],
                'formula': mod['formula'],
                'mass_avg': mod['mass_avg'],
                'reference_moiety': mod['reference_moiety'],
                'source': 'MODOMICS',
            }

            yield (f"MODOMICS:{mod['id']}", "RNAModification", props)
            count += 1

        logger.info(f"MODOMICS: Generated {count} RNAModification nodes")

    def get_edges(self):
        """No edges."""
        logger.info("MODOMICS: No edges to generate")
        return iter([])
