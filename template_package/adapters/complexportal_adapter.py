"""
ComplexPortal Adapter for BioCypher.

Loads curated protein complex data from ComplexPortal TSV files and generates:
- ProteinComplex nodes (macromolecular complexes)
- ComplexContainsProtein edges (complex → gene with stoichiometry)
"""

import re
import json
from pathlib import Path
from biocypher._logger import logger


class ComplexPortalAdapter:
    def __init__(self, data_dir="template_package/data/complexportal",
                 ortholog_file="template_package/mappings/mouse_to_human_orthologs.json"):
        self.data_dir = data_dir
        self.complexes = []
        self.orthologs = {}
        self._load_orthologs(ortholog_file)
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_orthologs(self, path):
        """Load mouse-to-human ortholog mappings."""
        try:
            with open(path, 'r') as f:
                self.orthologs = json.load(f)
            logger.info(f"ComplexPortal: Loaded {len(self.orthologs)} ortholog mappings")
        except FileNotFoundError:
            logger.warning(f"ComplexPortal: Ortholog file not found: {path}")

    def _parse_components(self, components_str):
        """
        Parse component string like 'P84022(1)|Q13485(1)|Q15796(1)'.
        Returns list of (uniprot_id, stoichiometry) tuples.
        """
        components = []
        if not components_str or components_str == '-':
            return components

        for part in components_str.split('|'):
            part = part.strip()
            if not part:
                continue
            # Match pattern: UNIPROT_ID(stoichiometry) or just UNIPROT_ID
            match = re.match(r'([A-Z0-9_-]+)\((\d+)\)', part)
            if match:
                uid = match.group(1)
                stoich = int(match.group(2))
                components.append((uid, stoich))
            else:
                # Just a bare ID
                components.append((part, 1))

        return components

    def _parse_go_annotations(self, go_str):
        """Parse GO annotation string like 'GO:0071144(name)|GO:0003690(name)'."""
        go_terms = []
        if not go_str or go_str == '-':
            return go_terms

        for part in go_str.split('|'):
            match = re.match(r'(GO:\d+)', part.strip())
            if match:
                go_terms.append(match.group(1))

        return go_terms

    def _load_data(self):
        """Load ComplexPortal TSV files for human and mouse."""
        for species_file, species_name, taxid in [
            ('homo_sapiens.tsv', 'Homo sapiens', '9606'),
            ('mus_musculus.tsv', 'Mus musculus', '10090'),
        ]:
            filepath = Path(self.data_dir) / species_file
            if not filepath.exists():
                logger.warning(f"ComplexPortal: {species_file} not found, skipping")
                continue

            logger.info(f"ComplexPortal: Loading {species_file}...")
            with open(filepath, 'r', encoding='utf-8') as f:
                header = None
                for line in f:
                    if line.startswith('#'):
                        # Parse header
                        header = line.lstrip('#').strip().split('\t')
                        continue
                    if not line.strip():
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    complex_data = {
                        'complex_ac': fields[0].strip(),
                        'name': self._sanitize(fields[1].strip()),
                        'aliases': self._sanitize(fields[2].strip()) if len(fields) > 2 else '',
                        'taxid': taxid,
                        'species': species_name,
                        'components_str': fields[4].strip() if len(fields) > 4 else '',
                        'go_annotations': self._parse_go_annotations(
                            fields[7].strip() if len(fields) > 7 else ''),
                        'description': self._sanitize(
                            fields[9].strip() if len(fields) > 9 else ''),
                        'complex_assembly': self._sanitize(
                            fields[11].strip() if len(fields) > 11 else ''),
                    }

                    # Parse components
                    complex_data['components'] = self._parse_components(
                        complex_data['components_str'])

                    self.complexes.append(complex_data)

        logger.info(f"ComplexPortal: Loaded {len(self.complexes)} complexes total")

    def _map_to_human(self, uniprot_id, species):
        """Map a protein ID to human ortholog if mouse."""
        if species == 'Mus musculus':
            return self.orthologs.get(uniprot_id, uniprot_id)
        return uniprot_id

    def get_nodes(self):
        """
        Generate ProteinComplex nodes.
        Yields: (id, label, properties)
        """
        logger.info("ComplexPortal: Generating nodes...")
        node_count = 0

        # Also collect unique Gene nodes from complex components
        gene_nodes = set()

        for cpx in self.complexes:
            props = {
                'name': cpx['name'],
                'aliases': cpx['aliases'],
                'description': cpx['description'],
                'complex_assembly': cpx['complex_assembly'],
                'species': cpx['species'],
                'go_annotations': cpx['go_annotations'],
                'num_components': len(cpx['components']),
            }

            yield (cpx['complex_ac'], "ProteinComplex", props)
            node_count += 1

            # Track gene nodes (don't yield them - they may already exist from LIANA)
            for uid, _ in cpx['components']:
                human_id = self._map_to_human(uid, cpx['species'])
                gene_nodes.add(human_id)

        logger.info(f"ComplexPortal: Generated {node_count} ProteinComplex nodes "
                     f"(referencing {len(gene_nodes)} unique genes)")

    def get_edges(self):
        """
        Generate ComplexContainsProtein edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("ComplexPortal: Generating edges...")
        edge_count = 0

        for cpx in self.complexes:
            for uid, stoich in cpx['components']:
                human_id = self._map_to_human(uid, cpx['species'])

                props = {
                    'stoichiometry': stoich,
                    'species': cpx['species'],
                }

                # If it was a mouse protein, track the original ID
                if cpx['species'] == 'Mus musculus' and human_id != uid:
                    props['original_id'] = uid

                yield (
                    None,
                    cpx['complex_ac'],
                    human_id,
                    "ComplexContainsProtein",
                    props
                )
                edge_count += 1

        logger.info(f"ComplexPortal: Generated {edge_count} ComplexContainsProtein edges")
