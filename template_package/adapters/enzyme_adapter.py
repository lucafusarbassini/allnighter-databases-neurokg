"""
ExPASy ENZYME Adapter for BioCypher.

Loads the ExPASy ENZYME nomenclature database and generates:
- Enzyme nodes (EC number classifications with catalytic activities)
- EnzymeCatalyzedBy edges (enzyme → protein via UniProt cross-references)

The ENZYME database contains nomenclature data for all characterized
enzymes, classified by EC numbers following IUBMB recommendations.
"""

from pathlib import Path
from biocypher._logger import logger


class EnzymeAdapter:
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
        """Load enzyme.dat file."""
        path = self.data_dir / 'enzyme.dat'
        if not path.exists():
            logger.warning("ENZYME: enzyme.dat not found")
            return

        logger.info("ENZYME: Loading enzyme nomenclature data...")
        count = 0
        current = {}

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('//'):
                    if current.get('id'):
                        self.enzymes.append(current)
                        count += 1
                    current = {}
                elif line.startswith('ID   '):
                    current['id'] = line[5:].strip()
                elif line.startswith('DE   '):
                    de = current.get('de', '')
                    current['de'] = (de + ' ' + line[5:].strip()).strip()
                elif line.startswith('AN   '):
                    ans = current.get('an', [])
                    ans.append(line[5:].strip().rstrip('.'))
                    current['an'] = ans
                elif line.startswith('CA   '):
                    cas = current.get('ca', [])
                    cas.append(line[5:].strip())
                    current['ca'] = cas
                elif line.startswith('DR   '):
                    # Parse UniProt cross-references
                    drs = current.get('dr', [])
                    # Format: "P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU; ..."
                    refs = line[5:].strip().rstrip(';').split(';')
                    for ref in refs:
                        ref = ref.strip()
                        if ',' in ref:
                            parts = ref.split(',')
                            acc = parts[0].strip()
                            entry = parts[1].strip().rstrip(';')
                            if '_HUMAN' in entry:
                                drs.append(acc)
                    current['dr'] = drs

        logger.info(f"ENZYME: Loaded {count} enzyme entries")

    def get_nodes(self):
        """
        Generate Enzyme nodes.
        Yields: (id, label, properties)
        """
        logger.info("ENZYME: Generating nodes...")
        count = 0

        for enzyme in self.enzymes:
            ec_id = enzyme['id']
            # Skip transferred/deleted entries
            if 'Transferred entry' in enzyme.get('de', '') or 'Deleted entry' in enzyme.get('de', ''):
                continue

            props = {
                'name': self._sanitize(enzyme.get('de', '')),
                'alternative_names': '|'.join(self._sanitize(a) for a in enzyme.get('an', [])[:5]),
                'catalytic_activity': self._sanitize(' '.join(enzyme.get('ca', []))[:500]),
                'ec_class': ec_id.split('.')[0] if ec_id else '',
                'source': 'ExPASy_ENZYME',
            }

            yield (f"EC:{ec_id}", "Enzyme", props)
            count += 1

        logger.info(f"ENZYME: Generated {count} Enzyme nodes")

    def get_edges(self):
        """
        Generate EnzymeCatalyzedBy edges (EC → UniProt human proteins).
        Yields: (id, source, target, label, properties)
        """
        logger.info("ENZYME: Generating edges...")
        count = 0

        for enzyme in self.enzymes:
            ec_id = enzyme['id']
            if 'Transferred entry' in enzyme.get('de', '') or 'Deleted entry' in enzyme.get('de', ''):
                continue

            for uniprot_acc in enzyme.get('dr', []):
                props = {
                    'ec_number': ec_id,
                    'source': 'ExPASy_ENZYME',
                }

                yield (
                    None,
                    uniprot_acc,
                    f"EC:{ec_id}",
                    "EnzymeCatalyzedBy",
                    props
                )
                count += 1

        logger.info(f"ENZYME: Generated {count} EnzymeCatalyzedBy edges")
