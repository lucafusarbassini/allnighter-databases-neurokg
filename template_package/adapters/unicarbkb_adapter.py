"""
UniCarbKB / GlyGen Glycan Adapter for BioCypher.

Loads glycan structure data from GlyGen (which integrates UniCarbKB data).
GlyGen is a comprehensive glycoinformatics resource that provides
standardized glycan data from GlyTouCan, UniCarbKB, and other repositories.

Generates:
- Glycan nodes (N-linked and other glycan structures)
"""

import csv
import json
from pathlib import Path
from biocypher._logger import logger


class UniCarbKBAdapter:
    def __init__(self, data_dir="template_package/data/unicarbkb"):
        self.data_dir = Path(data_dir)
        self.glycans = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load glycan data from GlyGen TSV files."""
        # Load N-linked glycans (most comprehensive)
        nlinked_path = self.data_dir / 'glygen_nlinked_glycans.tsv'
        if nlinked_path.exists():
            self._load_tsv(nlinked_path, 'N-linked')

        # Load general glycans
        general_path = self.data_dir / 'glygen_glycans.tsv'
        if general_path.exists():
            self._load_tsv(general_path, 'general')

        # Also try JSON files for richer data
        nlinked_json = self.data_dir / 'glygen_nlinked_glycans.json'
        if nlinked_json.exists():
            self._load_json(nlinked_json)

        logger.info(f"UniCarbKB/GlyGen: Loaded {len(self.glycans)} unique glycans")

    def _load_tsv(self, path, glycan_type):
        """Load glycan data from TSV file."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                count = 0
                for row in reader:
                    acc = row.get('glytoucan_ac', '').strip()
                    if not acc:
                        continue

                    if acc not in self.glycans:
                        self.glycans[acc] = {
                            'glytoucan_ac': acc,
                            'mass': row.get('mass', ''),
                            'glycan_type': glycan_type,
                            'byonic': row.get('byonic', ''),
                            'hit_score': row.get('hit_score', ''),
                            'publication_count': row.get('publication_count', '0'),
                        }
                        count += 1

            logger.info(f"UniCarbKB/GlyGen: Loaded {count} glycans from {path.name}")
        except Exception as e:
            logger.warning(f"UniCarbKB/GlyGen: Error loading {path.name}: {e}")

    def _load_json(self, path):
        """Load additional glycan data from JSON for richer properties."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            if isinstance(data, list):
                for item in data:
                    acc = item.get('glytoucan_ac', '')
                    if not acc:
                        continue
                    if acc in self.glycans:
                        # Enrich existing entry
                        if 'mass' in item and item['mass']:
                            self.glycans[acc]['mass'] = item['mass']
                    else:
                        self.glycans[acc] = {
                            'glytoucan_ac': acc,
                            'mass': item.get('mass', ''),
                            'glycan_type': 'N-linked',
                            'byonic': item.get('byonic', ''),
                            'hit_score': item.get('hit_score', ''),
                            'publication_count': str(item.get('publication_count', 0)),
                        }
        except Exception as e:
            logger.warning(f"UniCarbKB/GlyGen: Error loading JSON: {e}")

    def get_nodes(self):
        """
        Generate Glycan nodes.
        Yields: (id, label, properties)
        """
        logger.info("UniCarbKB/GlyGen: Generating glycan nodes...")
        count = 0

        for acc, glycan in self.glycans.items():
            mass_str = str(glycan.get('mass', '')).strip()
            try:
                mass = float(mass_str) if mass_str else 0.0
            except ValueError:
                mass = 0.0

            hit_score_str = str(glycan.get('hit_score', '')).strip()
            try:
                hit_score = float(hit_score_str) if hit_score_str else 0.0
            except ValueError:
                hit_score = 0.0

            pub_str = str(glycan.get('publication_count', '0')).strip()
            try:
                pub_count = int(pub_str) if pub_str else 0
            except ValueError:
                pub_count = 0

            props = {
                'composition': self._sanitize(glycan.get('byonic', '')),
                'mass': mass,
                'glycan_type': glycan.get('glycan_type', ''),
                'hit_score': hit_score,
                'publication_count': pub_count,
                'source': 'GlyGen/UniCarbKB',
            }

            yield (acc, "Glycan", props)
            count += 1

        logger.info(f"UniCarbKB/GlyGen: Generated {count} glycan nodes")

    def get_edges(self):
        """No edges - glycan-protein links from GlyGen require additional data."""
        logger.info("UniCarbKB/GlyGen: No edges to generate")
        return iter([])
