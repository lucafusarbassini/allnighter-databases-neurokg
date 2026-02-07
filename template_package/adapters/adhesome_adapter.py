"""
Adhesome Adapter for BioCypher.

Loads cell adhesion component and interaction data from Adhesome.org.
The adhesome is the collection of proteins involved in cell-matrix
and cell-cell adhesion signaling.

Note: Data files may contain download errors (e.g. "404: Not Found").
In that case, this adapter gracefully handles empty/invalid data.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class AdhesomeAdapter:
    def __init__(self, data_dir="template_package/data/adhesome"):
        self.data_dir = Path(data_dir)
        self.components = []
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
        """Load adhesome components and interactions."""
        # Load components
        comp_path = self.data_dir / 'components.csv'
        if comp_path.exists():
            try:
                with open(comp_path, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                    if content.startswith('404') or len(content) < 50:
                        logger.warning("Adhesome: components.csv appears invalid (likely download error)")
                    else:
                        f.seek(0)
                        reader = csv.DictReader(f)
                        for row in reader:
                            self.components.append(row)
                        logger.info(f"Adhesome: Loaded {len(self.components)} components")
            except Exception as e:
                logger.warning(f"Adhesome: Error reading components.csv: {e}")

        # Load interactions
        int_path = self.data_dir / 'interactions.csv'
        if int_path.exists():
            try:
                with open(int_path, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                    if content.startswith('404') or len(content) < 50:
                        logger.warning("Adhesome: interactions.csv appears invalid (likely download error)")
                    else:
                        f.seek(0)
                        reader = csv.DictReader(f)
                        for row in reader:
                            self.interactions.append(row)
                        logger.info(f"Adhesome: Loaded {len(self.interactions)} interactions")
            except Exception as e:
                logger.warning(f"Adhesome: Error reading interactions.csv: {e}")

        if not self.components and not self.interactions:
            logger.warning("Adhesome: No valid data loaded (files may have failed to download)")

    def get_nodes(self):
        """Generate adhesome component nodes (if data is valid)."""
        logger.info(f"Adhesome: Generating nodes from {len(self.components)} components...")
        count = 0

        for comp in self.components:
            gene_name = comp.get('Official Symbol', comp.get('Gene', '')).strip()
            if not gene_name:
                continue

            uniprot = comp.get('Swiss-Prot ID', comp.get('UniProt', '')).strip()
            node_id = uniprot if uniprot else gene_name

            props = {
                'gene_name': self._sanitize(gene_name),
                'functional_category': self._sanitize(comp.get('Functional Category', '')),
                'source': 'Adhesome',
            }

            yield (node_id, "Gene", props)
            count += 1

        logger.info(f"Adhesome: Generated {count} Gene nodes")

    def get_edges(self):
        """Generate adhesome interaction edges (if data is valid)."""
        logger.info(f"Adhesome: Generating edges from {len(self.interactions)} interactions...")
        count = 0

        for inter in self.interactions:
            source = inter.get('Source', '').strip()
            target = inter.get('Target', '').strip()
            if not source or not target:
                continue

            props = {
                'interaction_type': self._sanitize(inter.get('Type', '')),
                'source_db': 'Adhesome',
            }

            yield (None, source, target, "AdhesomeInteraction", props)
            count += 1

        logger.info(f"Adhesome: Generated {count} edges")
