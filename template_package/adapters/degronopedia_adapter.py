"""
Degronopedia Adapter for BioCypher.

Loads degron motif data and generates:
- Degron nodes (degradation signal motifs recognized by the UPS)

Degronopedia catalogs degron motifs - short protein sequences that
target proteins for ubiquitin-proteasome-dependent degradation.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class DegronopediaAdapter:
    def __init__(self, data_dir="template_package/data/degronopedia"):
        self.data_dir = Path(data_dir)
        self.degrons = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Degronopedia degron data."""
        path = self.data_dir / 'Degrons.tsv'
        if not path.exists():
            logger.warning("Degronopedia: TSV file not found")
            return

        logger.info("Degronopedia: Loading degron data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                degron = row.get('Degron', '').strip()
                regex = row.get('Degron_regex', '').strip()
                organism = row.get('Organism', '').strip()
                location = row.get('Degron_location', '').strip()
                degron_type = row.get('Degron_type', '').strip()
                ups_components = row.get('Known_UPS_components_recognizing_degron', '').strip()

                if not degron:
                    continue
                if 'sapiens' not in organism:
                    continue

                self.degrons.append({
                    'degron': degron,
                    'regex': regex,
                    'location': location,
                    'degron_type': degron_type,
                    'ups_components': ups_components,
                })
                count += 1

        logger.info(f"Degronopedia: Loaded {count} degron motifs")

    def get_nodes(self):
        """
        Generate Degron nodes.
        Yields: (id, label, properties)
        """
        logger.info("Degronopedia: Generating nodes...")
        count = 0
        seen = set()

        for deg in self.degrons:
            key = deg['degron']
            if key in seen:
                continue
            seen.add(key)

            props = {
                'regex': deg['regex'],
                'location': deg['location'],
                'degron_type': deg['degron_type'],
                'ups_components': self._sanitize(deg['ups_components']),
                'source': 'Degronopedia',
            }

            yield (f"DEGRON:{key}", "Degron", props)
            count += 1

        logger.info(f"Degronopedia: Generated {count} Degron nodes")

    def get_edges(self):
        """No edges."""
        logger.info("Degronopedia: No edges to generate")
        return iter([])
