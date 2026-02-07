"""
ArrestinDB Adapter for BioCypher.

Loads ArrestinDB GPCR-arrestin coupling data and generates:
- ArrestinCoupling edges (GPCR → arrestin coupling with selectivity)

ArrestinDB catalogs experimental data on GPCR-arrestin interactions,
providing quantitative coupling measurements.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class ArrestinDBAdapter:
    def __init__(self, data_dir="template_package/data/arrestindb"):
        self.data_dir = Path(data_dir)
        self.couplings = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ArrestinDB coupling data."""
        path = self.data_dir / 'arrestin_couplings_data.csv'
        if not path.exists():
            logger.warning("ArrestinDB: coupling data not found")
            return

        logger.info("ArrestinDB: Loading GPCR-arrestin coupling data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = None
            for row in reader:
                if header is None:
                    header = row
                    continue

                if len(row) < 5:
                    continue

                # Parse the CSV structure
                source = row[1] if len(row) > 1 else ''
                receptor = row[2] if len(row) > 2 else ''
                uniprot = row[4] if len(row) > 4 else ''
                ligand_name = row[6] if len(row) > 6 else ''

                if not receptor or not uniprot:
                    continue

                self.couplings.append({
                    'source_study': source.strip(),
                    'receptor': receptor.strip(),
                    'uniprot': uniprot.strip(),
                    'ligand': ligand_name.strip(),
                })
                count += 1

        logger.info(f"ArrestinDB: Loaded {count} coupling entries")

    def get_nodes(self):
        """No new nodes."""
        logger.info("ArrestinDB: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate ArrestinCoupling edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("ArrestinDB: Generating edges...")
        seen = set()
        count = 0

        for coupling in self.couplings:
            key = (coupling['uniprot'], coupling['receptor'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'receptor_name': self._sanitize(coupling['receptor']),
                'ligand': self._sanitize(coupling['ligand']),
                'source': 'ArrestinDB',
            }

            yield (
                None,
                coupling['uniprot'],
                "ARRESTIN_PATHWAY",
                "ArrestinCoupling",
                props
            )
            count += 1

        logger.info(f"ArrestinDB: Generated {count} ArrestinCoupling edges")
