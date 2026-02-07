"""
MetalPDB Adapter for BioCypher.

Loads MetalPDB metal binding site data and generates:
- MetalBindingSite edges (protein → metal binding with coordination details)

MetalPDB is a database of metal binding sites in 3D structures from PDB,
providing detailed coordination geometry information.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MetalPDBAdapter:
    def __init__(self, data_dir="template_package/data/metalpdb"):
        self.data_dir = Path(data_dir)
        self.sites = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load MetalPDB metal site data for human structures."""
        metals = ['Zn', 'Fe', 'Ca', 'Mg', 'Cu', 'Mn', 'Co', 'Ni']
        total = 0

        for metal in metals:
            path = self.data_dir / f'human_{metal}_sites.json'
            if not path.exists():
                continue

            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            count = 0
            for site in data:
                site_id = site.get('site', '').strip()
                pfam = site.get('pfam', '').strip()
                ec = site.get('ec_number', '')

                # Extract PDB ID from site name (e.g., "1au1_1" -> "1au1")
                pdb_id = site_id.split('_')[0] if site_id else ''

                # Get metal info
                metal_info = site.get('metals', [])
                num_ligands = sum(len(m.get('ligands', [])) for m in metal_info)

                self.sites.append({
                    'site_id': site_id,
                    'pdb_id': pdb_id,
                    'metal': metal,
                    'pfam': pfam,
                    'ec_number': ec or '',
                    'num_ligands': num_ligands,
                })
                count += 1

            total += count

        logger.info(f"MetalPDB: Loaded {total} metal binding sites across {len(metals)} metals")

    def get_nodes(self):
        """No new nodes."""
        logger.info("MetalPDB: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate MetalBindingSite edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("MetalPDB: Generating edges...")
        seen = set()
        count = 0

        for site in self.sites:
            # Deduplicate by site_id
            if site['site_id'] in seen:
                continue
            seen.add(site['site_id'])

            props = {
                'metal': site['metal'],
                'pfam': self._sanitize(site['pfam']),
                'ec_number': site['ec_number'],
                'num_ligands': site['num_ligands'],
                'source': 'MetalPDB',
            }

            yield (
                None,
                site['pdb_id'],
                f"METAL:{site['metal']}",
                "MetalBindingSite",
                props
            )
            count += 1

        logger.info(f"MetalPDB: Generated {count} MetalBindingSite edges")
