"""
SABIO-RK Adapter for BioCypher.

Loads SABIO-RK enzyme kinetic data and generates:
- KineticParameter edges (reaction kinetic parameters)

SABIO-RK is a database for biochemical reaction kinetics, providing
curated kinetic parameters (Km, Vmax, Ki, etc.) for biochemical reactions.
"""

import json
import xml.etree.ElementTree as ET
from pathlib import Path
from biocypher._logger import logger


class SABIORKAdapter:
    def __init__(self, data_dir="template_package/data/sabiork"):
        self.data_dir = Path(data_dir)
        self.kinetic_laws = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load SABIO-RK kinetic data from JSON and XML files."""
        # Load JSON kinetics
        json_path = self.data_dir / 'human_kinetics.json'
        if json_path.exists():
            logger.info("SABIO-RK: Loading JSON kinetic data...")
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            laws = data.get('kinetic_laws', [])
            for law in laws:
                reaction_id = law.get('reaction_id', '')
                reaction_name = law.get('reaction_name', '')
                params = law.get('parameters', [])

                for param in params:
                    self.kinetic_laws.append({
                        'reaction_id': reaction_id,
                        'reaction_name': reaction_name,
                        'param_name': param.get('name', ''),
                        'param_value': param.get('value', ''),
                        'param_units': param.get('units', ''),
                        'sbo_term': param.get('sboTerm', ''),
                    })

        # Also try to parse XML files for additional data
        for xml_file in sorted(self.data_dir.glob('human_kinetics_batch*.xml')):
            try:
                tree = ET.parse(xml_file)
                root = tree.getroot()
                ns = {'sbml': 'http://www.sbml.org/sbml/level2/version4'}
                for reaction in root.iter('{http://www.sbml.org/sbml/level2/version4}reaction'):
                    rxn_id = reaction.get('id', '')
                    rxn_name = reaction.get('name', '')
                    for param in reaction.iter('{http://www.sbml.org/sbml/level2/version4}parameter'):
                        self.kinetic_laws.append({
                            'reaction_id': rxn_id,
                            'reaction_name': rxn_name,
                            'param_name': param.get('name', ''),
                            'param_value': param.get('value', ''),
                            'param_units': param.get('units', ''),
                            'sbo_term': param.get('sboTerm', ''),
                        })
            except ET.ParseError:
                logger.warning(f"SABIO-RK: Could not parse {xml_file.name}")

        logger.info(f"SABIO-RK: Loaded {len(self.kinetic_laws)} kinetic parameters")

    def get_nodes(self):
        """No new nodes."""
        logger.info("SABIO-RK: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate KineticParameter edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("SABIO-RK: Generating edges...")
        count = 0

        for law in self.kinetic_laws:
            props = {
                'param_name': self._sanitize(law['param_name']),
                'param_value': law['param_value'],
                'param_units': law['param_units'],
                'sbo_term': law['sbo_term'],
                'source': 'SABIO-RK',
            }

            yield (
                None,
                law['reaction_id'],
                f"KINETIC:{law['reaction_id']}_{law['param_name']}",
                "KineticParameter",
                props
            )
            count += 1

        logger.info(f"SABIO-RK: Generated {count} KineticParameter edges")
