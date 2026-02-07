"""
BRENDA Enzyme Adapter for BioCypher.

Loads the ExPASy/BRENDA enzyme.dat flat file and generates:
- Enzyme nodes with detailed catalytic activity and classification
- Edges linking enzymes to proteins (via DR cross-references)

The enzyme.dat file follows the Swiss-Prot-like flat file format
with ID (EC number), DE (description), AN (alternative names),
CA (catalytic activity), CC (comments), DR (cross-references).
"""

import re
from pathlib import Path
from biocypher._logger import logger


class BRENDAAdapter:
    def __init__(self, data_dir="template_package/data/brenda"):
        self.data_dir = Path(data_dir)
        self.enzymes = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Parse the enzyme.dat flat file."""
        dat_path = self.data_dir / 'enzyme.dat'
        if not dat_path.exists():
            logger.warning("BRENDA: enzyme.dat file not found")
            return

        logger.info("BRENDA: Parsing enzyme.dat...")

        current = {
            'id': '', 'name': '', 'alt_names': [],
            'catalytic_activity': [], 'comments': [],
            'cross_refs': []
        }

        with open(dat_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')

                if line.startswith('//'):
                    # End of record
                    if current['id']:
                        self.enzymes.append({
                            'ec_number': current['id'],
                            'name': ' '.join(current.get('name_parts', [])).strip().rstrip('.'),
                            'alternative_names': '; '.join(current['alt_names']),
                            'catalytic_activity': ' '.join(current['catalytic_activity']).strip(),
                            'comments': ' '.join(current['comments']).strip(),
                            'cross_refs': current['cross_refs'],
                        })
                    current = {
                        'id': '', 'name_parts': [], 'alt_names': [],
                        'catalytic_activity': [], 'comments': [],
                        'cross_refs': []
                    }
                    continue

                if line.startswith('ID   '):
                    current['id'] = line[5:].strip()
                elif line.startswith('DE   '):
                    if 'name_parts' not in current:
                        current['name_parts'] = []
                    current['name_parts'].append(line[5:].strip())
                elif line.startswith('AN   '):
                    an = line[5:].strip().rstrip('.')
                    if an:
                        current['alt_names'].append(an)
                elif line.startswith('CA   '):
                    current['catalytic_activity'].append(line[5:].strip())
                elif line.startswith('CC   '):
                    cc_text = line[5:].strip()
                    if cc_text and not cc_text.startswith('---'):
                        current['comments'].append(cc_text)
                elif line.startswith('DR   '):
                    # DR lines contain cross-refs like: P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU;
                    dr_text = line[5:].strip()
                    refs = dr_text.split(';')
                    for ref in refs:
                        ref = ref.strip()
                        if ',' in ref:
                            parts = ref.split(',')
                            uniprot_id = parts[0].strip()
                            entry_name = parts[1].strip() if len(parts) > 1 else ''
                            if uniprot_id and '_HUMAN' in entry_name:
                                current['cross_refs'].append(uniprot_id)

        logger.info(f"BRENDA: Loaded {len(self.enzymes)} enzyme entries")

    def get_nodes(self):
        """
        Generate Enzyme nodes.
        Yields: (id, label, properties)
        """
        logger.info("BRENDA: Generating enzyme nodes...")
        count = 0

        for enz in self.enzymes:
            ec = enz['ec_number']
            if not ec or ec.startswith('Transferred') or ec.startswith('Deleted'):
                continue

            # Determine EC class
            ec_parts = ec.split('.')
            ec_class = ''
            class_map = {
                '1': 'Oxidoreductases', '2': 'Transferases',
                '3': 'Hydrolases', '4': 'Lyases',
                '5': 'Isomerases', '6': 'Ligases',
                '7': 'Translocases'
            }
            if ec_parts and ec_parts[0] in class_map:
                ec_class = class_map[ec_parts[0]]

            props = {
                'name': self._sanitize(enz['name']),
                'alternative_names': self._sanitize(enz['alternative_names']),
                'catalytic_activity': self._sanitize(enz['catalytic_activity']),
                'ec_class': ec_class,
                'source': 'BRENDA/ExPASy',
            }

            yield (f"EC:{ec}", "Enzyme", props)
            count += 1

        logger.info(f"BRENDA: Generated {count} Enzyme nodes")

    def get_edges(self):
        """
        Generate edges linking enzymes to human proteins.
        Yields: (id, source, target, label, properties)
        """
        logger.info("BRENDA: Generating enzyme-protein edges...")
        count = 0

        for enz in self.enzymes:
            ec = enz['ec_number']
            if not ec or ec.startswith('Transferred') or ec.startswith('Deleted'):
                continue

            for uniprot_id in enz['cross_refs']:
                props = {
                    'ec_number': ec,
                    'source': 'BRENDA/ExPASy',
                }
                yield (None, f"EC:{ec}", uniprot_id, "EnzymeCatalyzedBy", props)
                count += 1

        logger.info(f"BRENDA: Generated {count} enzyme-protein edges")
