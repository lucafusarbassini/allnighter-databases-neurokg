"""
Membranome Adapter for BioCypher.

Loads Membranome single-pass transmembrane protein data and generates:
- MembraneProtein nodes (transmembrane proteins with topology)

Membranome catalogs single-pass transmembrane proteins with their
structural and topological properties.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MembranomeAdapter:
    def __init__(self, data_dir="template_package/data/membranome"):
        self.data_dir = Path(data_dir)
        self.proteins = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Membranome human protein data."""
        path = self.data_dir / 'human_proteins_all.json'
        if not path.exists():
            logger.warning("Membranome: data file not found")
            return

        logger.info("Membranome: Loading transmembrane protein data...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        objects = data.get('objects', [])
        count = 0
        seen = set()

        for obj in objects:
            prot_id = obj.get('id', '')
            name = obj.get('name', '').strip()
            uniprot_id = obj.get('uniprot_id', '').strip()
            thickness = obj.get('thickness', 0)
            tilt = obj.get('tilt', 0)

            if not name:
                continue

            # Deduplicate by uniprot_id (data contains pagination duplicates)
            dedup_key = uniprot_id or name
            if dedup_key in seen:
                continue
            seen.add(dedup_key)

            self.proteins.append({
                'id': str(prot_id),
                'name': name,
                'uniprot_entry': uniprot_id,
                'thickness': thickness or 0,
                'tilt': tilt or 0,
            })
            count += 1

        logger.info(f"Membranome: Loaded {count} transmembrane proteins (deduplicated)")

    def get_nodes(self):
        """
        Generate MembraneProtein nodes.
        Yields: (id, label, properties)
        """
        logger.info("Membranome: Generating nodes...")
        count = 0

        for prot in self.proteins:
            props = {
                'name': self._sanitize(prot['name']),
                'uniprot_entry': prot['uniprot_entry'],
                'membrane_thickness': prot['thickness'],
                'tilt_angle': prot['tilt'],
                'source': 'Membranome',
            }

            yield (f"MEMBRANOME:{prot['id']}", "MembraneProtein", props)
            count += 1

        logger.info(f"Membranome: Generated {count} MembraneProtein nodes")

    def get_edges(self):
        """No edges for now."""
        logger.info("Membranome: No edges to generate")
        return iter([])
