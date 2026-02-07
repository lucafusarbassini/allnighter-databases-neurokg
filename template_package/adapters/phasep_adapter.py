"""
PhaSepDB/PhaseP Adapter for BioCypher.

Loads phase separation protein data and generates:
- PhaseSeparatingProtein edges (protein → phase separation behavior)

PhaSepDB catalogs proteins that undergo liquid-liquid phase separation,
forming membraneless organelles and biomolecular condensates.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class PhasePAdapter:
    def __init__(self, data_dir="template_package/data/phasep"):
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
        """Load PhaseP human phase separation data."""
        path = self.data_dir / 'human_phase_separation_proteins.tsv'
        if not path.exists():
            logger.warning("PhaseP: TSV file not found")
            return

        logger.info("PhaseP: Loading phase separation data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                uniprot_id = row.get('uniprot_id', '').strip()
                organism = row.get('organism', '').strip()
                location = row.get('location', '').strip()
                material_state = row.get('material_state', '').strip()
                ps_class = row.get('class_', '').strip()
                mlo = row.get('mlo_normalized', '').strip()

                if not uniprot_id:
                    continue
                if 'sapiens' not in organism:
                    continue

                self.proteins.append({
                    'uniprot_id': uniprot_id,
                    'location': location,
                    'material_state': material_state,
                    'ps_class': ps_class,
                    'mlo': mlo,
                })
                count += 1

        logger.info(f"PhaseP: Loaded {count} human phase-separating proteins")

    def get_nodes(self):
        """No new nodes."""
        logger.info("PhaseP: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate PhaseSeparation edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PhaseP: Generating edges...")
        seen = set()
        count = 0

        for prot in self.proteins:
            if prot['uniprot_id'] in seen:
                continue
            seen.add(prot['uniprot_id'])

            props = {
                'location': self._sanitize(prot['location']),
                'material_state': prot['material_state'],
                'ps_class': prot['ps_class'],
                'mlo': self._sanitize(prot['mlo']),
                'source': 'PhaSepDB',
            }

            yield (
                None,
                prot['uniprot_id'],
                f"CONDENSATE:{prot['location'] or 'unknown'}",
                "PhaseSeparation",
                props
            )
            count += 1

        logger.info(f"PhaseP: Generated {count} PhaseSeparation edges")
