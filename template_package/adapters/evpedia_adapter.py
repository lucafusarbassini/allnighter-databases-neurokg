"""
EVpedia Adapter for BioCypher.

Loads extracellular vesicle (EV) proteomics and RNA data from EVpedia.
EVpedia is a community web portal for extracellular vesicle research.

Generates:
- EV origin nodes (tissue/cell sources of EVs)
- circRNA annotation edges (circRNAs detected in EVs)
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class EVpediaAdapter:
    def __init__(self, data_dir="template_package/data/evpedia"):
        self.data_dir = Path(data_dir)
        self.origins = []
        self.circrnas = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load EVpedia data files."""
        # Load EV origin/tissue data
        origin_path = self.data_dir / 'browse_origin.csv'
        if origin_path.exists():
            try:
                with open(origin_path, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        name = row.get('Tissue/Cell name', '').strip()
                        if not name:
                            continue
                        self.origins.append(row)
                logger.info(f"EVpedia: Loaded {len(self.origins)} EV tissue/cell origins")
            except Exception as e:
                logger.warning(f"EVpedia: Error reading browse_origin.csv: {e}")

        # Load circRNA annotations
        circ_path = self.data_dir / 'circRNAs_anno.csv'
        if circ_path.exists():
            try:
                count = 0
                with open(circ_path, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        circ_id = row.get('circID', '').strip().strip('"')
                        if not circ_id:
                            continue
                        self.circrnas.append({
                            'circ_id': circ_id,
                            'circbase_id': row.get('circBase ID', '').strip().strip('"'),
                            'position': row.get('Genomic position', '').strip().strip('"'),
                            'strand': row.get('Strand', '').strip().strip('"'),
                            'gene_symbol': row.get('Gene symbol', '').strip().strip('"'),
                            'gene_type': row.get('Gene type', '').strip().strip('"'),
                            'sample_type': row.get('Sample type', '').strip().strip('"'),
                        })
                        count += 1
                logger.info(f"EVpedia: Loaded {count} circRNA annotations")
            except Exception as e:
                logger.warning(f"EVpedia: Error reading circRNAs_anno.csv: {e}")

    def get_nodes(self):
        """
        Generate EV origin tissue/cell nodes.
        Yields: (id, label, properties)
        """
        logger.info("EVpedia: Generating EV origin nodes...")
        count = 0

        for origin in self.origins:
            name = origin.get('Tissue/Cell name', '').strip()
            full_name = origin.get('Full name', name).strip()
            category = origin.get('Main category', '').strip()
            cell_type = origin.get('Tissue/Cell type', '').strip()

            node_id = f"EVP:{name}"

            props = {
                'name': self._sanitize(full_name),
                'category': self._sanitize(category),
                'cell_type': self._sanitize(cell_type),
                'source': 'EVpedia',
            }

            yield (node_id, "EVOrigin", props)
            count += 1

        # Also create circRNA nodes
        seen_circ = set()
        for circ in self.circrnas:
            cid = circ['circ_id']
            if cid in seen_circ:
                continue
            seen_circ.add(cid)

            props = {
                'position': self._sanitize(circ['position']),
                'strand': circ['strand'],
                'gene_symbol': self._sanitize(circ['gene_symbol']),
                'gene_type': self._sanitize(circ['gene_type']),
                'sample_type': self._sanitize(circ['sample_type']),
                'source': 'EVpedia',
            }

            yield (f"EVP:{cid}", "CircularRNA", props)
            count += 1

        logger.info(f"EVpedia: Generated {count} nodes")

    def get_edges(self):
        """No edges generated currently."""
        logger.info("EVpedia: No edges to generate")
        return iter([])
