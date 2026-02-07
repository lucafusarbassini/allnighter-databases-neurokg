"""
Dfam Adapter for BioCypher.

Loads Dfam transposable element (TE) family data and generates:
- TransposableElementFamily nodes (TE families with classification)

Dfam is a database of repetitive DNA elements including transposable
elements, satellites, and other repeat families.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class DfamAdapter:
    def __init__(self, data_dir="template_package/data/dfam"):
        self.data_dir = Path(data_dir)
        self.families = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Dfam human TE families TSV."""
        path = self.data_dir / 'dfam_human_te_families.tsv'
        if not path.exists():
            logger.warning("Dfam: TE families file not found")
            return

        logger.info("Dfam: Loading transposable element families...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                accession = row.get('accession', '').strip()
                name = row.get('name', '').strip()
                title = row.get('title', '').strip()
                length = row.get('length', '0').strip()
                repeat_type = row.get('repeat_type', '').strip()
                repeat_subtype = row.get('repeat_subtype', '').strip()
                classification = row.get('classification', '').strip()

                if not accession:
                    continue

                self.families.append({
                    'accession': accession,
                    'name': name,
                    'title': title,
                    'length': int(length) if length.isdigit() else 0,
                    'repeat_type': repeat_type,
                    'repeat_subtype': repeat_subtype,
                    'classification': classification,
                })
                count += 1

        logger.info(f"Dfam: Loaded {count} TE families")

    def get_nodes(self):
        """
        Generate TransposableElementFamily nodes.
        Yields: (id, label, properties)
        """
        logger.info("Dfam: Generating nodes...")
        count = 0

        for fam in self.families:
            props = {
                'name': self._sanitize(fam['name']),
                'title': self._sanitize(fam['title']),
                'consensus_length': fam['length'],
                'repeat_type': fam['repeat_type'],
                'repeat_subtype': fam['repeat_subtype'],
                'classification': self._sanitize(fam['classification']),
                'source': 'Dfam',
            }

            yield (fam['accession'], "TransposableElementFamily", props)
            count += 1

        logger.info(f"Dfam: Generated {count} TransposableElementFamily nodes")

    def get_edges(self):
        """No edges for now."""
        logger.info("Dfam: No edges to generate")
        return iter([])
