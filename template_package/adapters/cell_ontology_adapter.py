"""
Cell Ontology (CL) Adapter for BioCypher.

Parses the Cell Ontology OBO file and generates:
- CellType nodes (cell type terms with definitions)
- CellTypeIsA edges (parent-child hierarchy)
"""

import re
from pathlib import Path
from biocypher._logger import logger


class CellOntologyAdapter:
    def __init__(self, data_dir="template_package/data/cell_ontology"):
        self.data_dir = Path(data_dir)
        self.terms = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = re.sub(r'<[^>]+>', '', text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load OBO file."""
        # Try files in order of preference
        for fname in ['cl-basic-purl.obo', 'cl-basic.obo', 'cl.obo', 'cl-purl.obo']:
            obo_path = self.data_dir / fname
            if obo_path.exists() and obo_path.stat().st_size > 100:
                logger.info(f"CellOntology: Loading {fname}...")
                self._parse_obo(obo_path)
                if self.terms:
                    break

        logger.info(f"CellOntology: Loaded {len(self.terms)} cell type terms")

    def _parse_obo(self, path):
        """Parse OBO format file."""
        current_term = None
        in_term = False

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip()

                if line == '[Term]':
                    if current_term and current_term.get('id'):
                        self.terms[current_term['id']] = current_term
                    current_term = {
                        'id': '',
                        'name': '',
                        'definition': '',
                        'synonyms': [],
                        'xrefs': [],
                        'is_a': [],
                        'is_obsolete': False,
                    }
                    in_term = True
                    continue
                elif line.startswith('['):
                    # Another stanza type (e.g., [Typedef])
                    if current_term and current_term.get('id'):
                        self.terms[current_term['id']] = current_term
                    current_term = None
                    in_term = False
                    continue

                if not in_term or current_term is None:
                    continue

                if line.startswith('id: '):
                    current_term['id'] = line[4:].strip()
                elif line.startswith('name: '):
                    current_term['name'] = line[6:].strip()
                elif line.startswith('def: '):
                    # Parse definition: "text" [refs]
                    match = re.match(r'def:\s+"(.*?)"\s*\[', line)
                    if match:
                        current_term['definition'] = match.group(1)
                elif line.startswith('synonym: '):
                    match = re.match(r'synonym:\s+"(.*?)"\s+', line)
                    if match:
                        current_term['synonyms'].append(match.group(1))
                elif line.startswith('xref: '):
                    xref = line[6:].strip().split(' ')[0]
                    current_term['xrefs'].append(xref)
                elif line.startswith('is_a: '):
                    parent = line[6:].strip().split(' ')[0]
                    current_term['is_a'].append(parent)
                elif line.startswith('is_obsolete: true'):
                    current_term['is_obsolete'] = True

        # Don't forget the last term
        if current_term and current_term.get('id'):
            self.terms[current_term['id']] = current_term

    def get_nodes(self):
        """
        Generate cell type nodes.
        Yields: (id, label, properties)
        """
        logger.info("CellOntology: Generating nodes...")
        count = 0

        for term_id, term in self.terms.items():
            if term.get('is_obsolete'):
                continue
            if not term_id.startswith('CL:'):
                continue

            props = {
                'name': self._sanitize(term.get('name', '')),
                'definition': self._sanitize(term.get('definition', '')),
                'synonyms': '|'.join(self._sanitize(s) for s in term.get('synonyms', [])),
                'xrefs': '|'.join(term.get('xrefs', [])),
                'source': 'CellOntology',
            }

            yield (term_id, "CellType", props)
            count += 1

        logger.info(f"CellOntology: Generated {count} CellType nodes")

    def get_edges(self):
        """
        Generate cell type hierarchy edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("CellOntology: Generating edges...")
        count = 0

        for term_id, term in self.terms.items():
            if term.get('is_obsolete'):
                continue
            if not term_id.startswith('CL:'):
                continue

            for parent_id in term.get('is_a', []):
                if parent_id.startswith('CL:'):
                    yield (None, term_id, parent_id, "CellTypeIsA", {})
                    count += 1

        logger.info(f"CellOntology: Generated {count} CellTypeIsA edges")
