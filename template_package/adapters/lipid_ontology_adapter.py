"""
LION (Lipid Ontology) Adapter for BioCypher.

Loads LION ontology terms and lipid-term associations from the
Lipid Ontology project (http://www.lipidontology.com/).

LION provides a controlled vocabulary for lipid biological functions,
physical properties, and classifications, linking individual lipid
species to biological concepts.

Generates:
- LIONTerm nodes (ontology terms describing lipid properties/functions)
- Lipid-to-LION-term association edges
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class LipidOntologyAdapter:
    def __init__(self, data_dir="template_package/data/lipid_ontology"):
        self.data_dir = Path(data_dir)
        self.terms = {}
        self.associations = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load LION terms and lipid-term associations."""
        # Load LION terms
        terms_path = self.data_dir / 'LION-terms.csv'
        if terms_path.exists():
            try:
                with open(terms_path, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        name = row.get('name', '').strip()
                        lion_id = row.get('LION', '').strip()
                        if lion_id:
                            self.terms[lion_id] = name
                logger.info(f"LION: Loaded {len(self.terms)} ontology terms")
            except Exception as e:
                logger.warning(f"LION: Error reading LION-terms.csv: {e}")

        # Load lipid-LION associations (wide matrix format)
        assoc_path = self.data_dir / 'all-LION-lipid-associations.csv'
        if assoc_path.exists():
            try:
                with open(assoc_path, 'r', encoding='utf-8') as f:
                    reader = csv.reader(f)
                    header = next(reader)  # First row: empty + LION IDs
                    lion_ids = header[1:]  # Skip first empty column

                    # Second row: empty + LION term names
                    names_row = next(reader)

                    count = 0
                    for row in reader:
                        if not row or not row[0]:
                            continue
                        lipid_name = row[0].strip()

                        for i, val in enumerate(row[1:], 0):
                            val = val.strip()
                            if val == 'x' and i < len(lion_ids):
                                lion_id = lion_ids[i].strip()
                                self.associations.append({
                                    'lipid_name': lipid_name,
                                    'lion_id': lion_id,
                                })
                                count += 1

                logger.info(f"LION: Loaded {count} lipid-term associations")
            except Exception as e:
                logger.warning(f"LION: Error reading associations: {e}")

    def get_nodes(self):
        """
        Generate LION term nodes.
        Yields: (id, label, properties)
        """
        logger.info("LION: Generating ontology term nodes...")
        count = 0

        for lion_id, name in self.terms.items():
            props = {
                'name': self._sanitize(name),
                'source': 'LION',
            }

            yield (lion_id, "LipidOntologyTerm", props)
            count += 1

        logger.info(f"LION: Generated {count} ontology term nodes")

    def get_edges(self):
        """
        Generate lipid-to-LION-term association edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info(f"LION: Generating edges from {len(self.associations)} associations...")
        count = 0

        for assoc in self.associations:
            lipid_name = assoc['lipid_name']
            lion_id = assoc['lion_id']

            # Use lipid name as source ID (will need to cross-ref with LIPIDMAPS)
            lipid_id = f"LION_LIPID:{lipid_name}"

            props = {
                'lipid_name': self._sanitize(lipid_name),
                'source': 'LION',
            }

            yield (None, lipid_id, lion_id, "LipidHasOntologyTerm", props)
            count += 1

        logger.info(f"LION: Generated {count} lipid-term edges")
