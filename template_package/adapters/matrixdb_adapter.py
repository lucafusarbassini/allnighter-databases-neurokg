"""
MatrixDB Adapter for BioCypher.

Loads MatrixDB extracellular matrix interaction data and generates:
- ECMInteraction edges (extracellular matrix protein-protein interactions)

MatrixDB is a database focused on interactions established by extracellular
matrix proteins, proteoglycans, and polysaccharides.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class MatrixDBAdapter:
    def __init__(self, data_dir="template_package/data/matrixdb"):
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
        """Load MatrixDB interaction data from MITAB format."""
        path = self.data_dir / 'matrixdb_CORE.tab'
        if not path.exists():
            logger.warning("MatrixDB: CORE interaction file not found")
            return

        logger.info("MatrixDB: Loading ECM interactions...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 15:
                    continue

                id_a = parts[0].strip()
                id_b = parts[1].strip()
                detection_method = parts[6].strip()
                publication = parts[8].strip()
                taxid_a = parts[9].strip()
                taxid_b = parts[10].strip()
                interaction_type = parts[11].strip()

                # Extract UniProt IDs
                a_id = id_a.split(':')[-1] if ':' in id_a else id_a
                b_id = id_b.split(':')[-1] if ':' in id_b else id_b

                # Clean up IDs (remove quotes)
                a_id = a_id.strip('"').strip("'")
                b_id = b_id.strip('"').strip("'")

                if not a_id or not b_id:
                    continue

                # Extract short detection method
                method_short = ''
                if 'MI:' in detection_method:
                    import re
                    m = re.search(r'"([^"]*)"[)]\s*$', detection_method)
                    if m:
                        method_short = m.group(1)

                self.interactions.append({
                    'id_a': a_id,
                    'id_b': b_id,
                    'detection_method': method_short,
                    'interaction_type': interaction_type,
                })
                count += 1

        logger.info(f"MatrixDB: Loaded {count} ECM interactions")

    def get_nodes(self):
        """No new nodes."""
        logger.info("MatrixDB: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate ECMInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("MatrixDB: Generating edges...")
        seen = set()
        count = 0

        for interaction in self.interactions:
            key = tuple(sorted([interaction['id_a'], interaction['id_b']]))
            if key in seen:
                continue
            seen.add(key)

            props = {
                'detection_method': self._sanitize(interaction['detection_method']),
                'source': 'MatrixDB',
            }

            yield (
                None,
                interaction['id_a'],
                interaction['id_b'],
                "ECMInteraction",
                props
            )
            count += 1

        logger.info(f"MatrixDB: Generated {count} ECMInteraction edges")
