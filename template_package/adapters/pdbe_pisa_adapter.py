"""
PDBe PISA (Protein Interfaces, Surfaces and Assemblies) Adapter for BioCypher.

Loads PDBe PISA interface data and generates:
- ProteinInterface edges (protein-protein interface with structural details)

PISA analyzes macromolecular interfaces in crystal structures from PDB,
identifying biologically relevant assemblies and interfaces.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from biocypher._logger import logger


class PDBePISAAdapter:
    def __init__(self, data_dir="template_package/data/pdbe_pisa"):
        self.data_dir = Path(data_dir)
        self.interfaces = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PDBe PISA interface XML files."""
        logger.info("PDBe PISA: Loading interface data...")
        count = 0

        for xml_file in sorted(self.data_dir.glob('pisa_interfaces_*.xml')):
            try:
                tree = ET.parse(xml_file)
                root = tree.getroot()

                pdb_entry = root.find('pdb_entry')
                if pdb_entry is None:
                    continue

                pdb_code = ''
                pdb_code_el = pdb_entry.find('pdb_code')
                if pdb_code_el is not None and pdb_code_el.text:
                    pdb_code = pdb_code_el.text.strip()

                # Process interfaces
                for interface in pdb_entry.iter('interface'):
                    iface_id = interface.findtext('id', '').strip()
                    area = interface.findtext('int_area', '0').strip()
                    solv_en = interface.findtext('int_solv_en', '0').strip()
                    n_hbonds = interface.findtext('n_h-bonds', '0').strip()
                    n_salt = interface.findtext('n_salt-bridges', '0').strip()

                    # Get molecule info
                    molecules = interface.findall('molecule')
                    chain_ids = []
                    for mol in molecules[:2]:
                        chain_id = mol.findtext('chain_id', '').strip()
                        if chain_id:
                            chain_ids.append(chain_id)

                    if len(chain_ids) >= 2:
                        try:
                            area_val = float(area)
                        except ValueError:
                            area_val = 0.0

                        try:
                            hbonds = int(float(n_hbonds))
                        except ValueError:
                            hbonds = 0

                        try:
                            salt_bridges = int(float(n_salt))
                        except ValueError:
                            salt_bridges = 0

                        self.interfaces.append({
                            'pdb_code': pdb_code,
                            'interface_id': iface_id,
                            'chain_a': chain_ids[0],
                            'chain_b': chain_ids[1],
                            'interface_area': area_val,
                            'n_hbonds': hbonds,
                            'n_salt_bridges': salt_bridges,
                        })
                        count += 1

            except ET.ParseError:
                logger.warning(f"PDBe PISA: Could not parse {xml_file.name}")

        logger.info(f"PDBe PISA: Loaded {count} protein interfaces from {len(list(self.data_dir.glob('pisa_interfaces_*.xml')))} structures")

    def get_nodes(self):
        """No new nodes."""
        logger.info("PDBe PISA: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate ProteinInterface edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PDBe PISA: Generating edges...")
        count = 0

        for iface in self.interfaces:
            source_id = f"{iface['pdb_code']}_{iface['chain_a']}"
            target_id = f"{iface['pdb_code']}_{iface['chain_b']}"

            props = {
                'pdb_code': iface['pdb_code'],
                'interface_area': iface['interface_area'],
                'n_hbonds': iface['n_hbonds'],
                'n_salt_bridges': iface['n_salt_bridges'],
                'source': 'PDBe_PISA',
            }

            yield (
                None,
                source_id,
                target_id,
                "ProteinInterface",
                props
            )
            count += 1

        logger.info(f"PDBe PISA: Generated {count} ProteinInterface edges")
