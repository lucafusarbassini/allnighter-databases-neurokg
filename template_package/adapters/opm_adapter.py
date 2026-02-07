"""
OPM Adapter for BioCypher.

Loads OPM (Orientations of Proteins in Membranes) data and generates:
- MembraneOrientation edges (protein → membrane orientation with tilt/thickness)

OPM provides spatial positions of membrane proteins relative to the lipid
bilayer, including tilt angle and membrane thickness.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class OPMAdapter:
    def __init__(self, data_dir="template_package/data/opm"):
        self.data_dir = Path(data_dir)
        self.structures = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load OPM membrane protein orientation data."""
        logger.info("OPM: Loading membrane protein orientations...")
        count = 0

        for fname in sorted(self.data_dir.glob('opm_structures*.json')):
            try:
                with open(fname, 'r', encoding='utf-8') as f:
                    data = json.load(f)

                objects = data.get('objects', [])
                for obj in objects:
                    pdbid = obj.get('pdbid', '')
                    if not pdbid:
                        continue

                    self.structures.append({
                        'pdbid': pdbid,
                        'name': obj.get('name', ''),
                        'thickness': obj.get('thickness', 0),
                        'tilt': obj.get('tilt', 0),
                        'gibbs': obj.get('gibbs', 0),
                        'resolution': obj.get('resolution', ''),
                        'subunit_segments': obj.get('subunit_segments', 0),
                        'family_name': obj.get('family_name_cache', ''),
                        'species_name': obj.get('species_name_cache', ''),
                        'membrane_name': obj.get('membrane_name_cache', ''),
                    })
                    count += 1
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"OPM: Error loading {fname.name}: {e}")

        logger.info(f"OPM: Loaded {count} membrane protein structures")

    def get_nodes(self):
        """No new nodes."""
        logger.info("OPM: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate MembraneOrientation edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("OPM: Generating edges...")
        seen = set()
        count = 0

        for struct in self.structures:
            pdbid = struct['pdbid']
            if pdbid in seen:
                continue
            seen.add(pdbid)

            props = {
                'name': self._sanitize(struct['name']),
                'thickness': struct['thickness'],
                'tilt': struct['tilt'],
                'gibbs_energy': struct['gibbs'],
                'resolution': struct['resolution'],
                'subunit_segments': struct['subunit_segments'],
                'family_name': self._sanitize(struct['family_name']),
                'species': self._sanitize(struct['species_name']),
                'membrane_type': self._sanitize(struct['membrane_name']),
                'source': 'OPM',
            }

            yield (
                None,
                f"PDB:{pdbid}",
                struct['membrane_name'],
                "MembraneOrientation",
                props
            )
            count += 1

        logger.info(f"OPM: Generated {count} MembraneOrientation edges")
