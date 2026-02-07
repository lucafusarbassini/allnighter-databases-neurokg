"""
InterPro Adapter for BioCypher.

Loads InterPro protein family/domain data and generates:
- InterProEntry nodes (protein families, domains, sites)
- InterProHierarchy edges (parent-child relationships between entries)
- InterProToGO edges (InterPro entries → GO terms)
"""

import re
from pathlib import Path
from biocypher._logger import logger


class InterProAdapter:
    def __init__(self, data_dir="template_package/data/interpro"):
        self.data_dir = data_dir
        self.entries = {}  # IPR ID -> {type, name, short_name, go_terms: []}
        self.hierarchy = []  # list of (parent, child)
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load all InterPro data files."""
        self._load_entry_list()
        self._load_names()
        self._load_hierarchy()
        self._load_interpro2go()

    def _load_entry_list(self):
        """Load entry.list file with entry types."""
        filepath = Path(self.data_dir) / 'entry.list'
        if not filepath.exists():
            logger.warning("InterPro: entry.list not found")
            return

        logger.info("InterPro: Loading entry list...")
        with open(filepath, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    ipr_id = parts[0].strip()
                    entry_type = parts[1].strip()
                    entry_name = self._sanitize(parts[2].strip())
                    self.entries[ipr_id] = {
                        'type': entry_type,
                        'name': entry_name,
                        'short_name': '',
                        'go_terms': [],
                    }

        logger.info(f"InterPro: Loaded {len(self.entries)} entries")

    def _load_names(self):
        """Load short names from names.dat."""
        filepath = Path(self.data_dir) / 'names.dat'
        if not filepath.exists():
            return

        logger.info("InterPro: Loading short names...")
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    ipr_id = parts[0].strip()
                    short_name = self._sanitize(parts[1].strip())
                    if ipr_id in self.entries:
                        self.entries[ipr_id]['short_name'] = short_name

    def _load_hierarchy(self):
        """Load parent-child hierarchy from ParentChildTreeFile.txt."""
        filepath = Path(self.data_dir) / 'ParentChildTreeFile.txt'
        if not filepath.exists():
            logger.warning("InterPro: ParentChildTreeFile.txt not found")
            return

        logger.info("InterPro: Loading hierarchy...")
        # The file uses '--' indentation to show hierarchy
        # IPR000008::C2 domain::
        # --IPR014705::Synaptotagmin-17, C2B domain::
        # ----IPR030537::...
        parent_stack = []  # stack of (indent_level, ipr_id)

        with open(filepath, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue

                # Count leading '--' indentation
                indent = 0
                while line[indent * 2:indent * 2 + 2] == '--':
                    indent += 1

                # Extract IPR ID
                content = line[indent * 2:]
                match = re.match(r'(IPR\d+)', content)
                if not match:
                    continue

                ipr_id = match.group(1)

                # Find parent at the right level
                while parent_stack and parent_stack[-1][0] >= indent:
                    parent_stack.pop()

                if parent_stack:
                    parent_id = parent_stack[-1][1]
                    self.hierarchy.append((parent_id, ipr_id))

                parent_stack.append((indent, ipr_id))

        logger.info(f"InterPro: Loaded {len(self.hierarchy)} hierarchy relationships")

    def _load_interpro2go(self):
        """Load InterPro to GO term mappings."""
        filepath = Path(self.data_dir) / 'interpro2go'
        if not filepath.exists():
            return

        logger.info("InterPro: Loading InterPro2GO mappings...")
        count = 0
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('!') or not line.strip():
                    continue
                # Format: InterPro:IPR000001 Kringle > GO:molecular_function ; GO:0003674
                match = re.match(r'InterPro:(IPR\d+).*>(.*?;\s*(GO:\d+))', line)
                if match:
                    ipr_id = match.group(1)
                    go_id = match.group(3)
                    if ipr_id in self.entries:
                        self.entries[ipr_id]['go_terms'].append(go_id)
                        count += 1

        logger.info(f"InterPro: Loaded {count} GO mappings")

    def get_nodes(self):
        """
        Generate InterProEntry nodes.
        Yields: (id, label, properties)
        """
        logger.info("InterPro: Generating nodes...")
        node_count = 0

        # Map InterPro types to BioLink-compatible labels
        type_map = {
            'Family': 'ProteinFamily',
            'Domain': 'ProteinDomain',
            'Repeat': 'ProteinDomain',
            'Homologous_superfamily': 'ProteinFamily',
            'Active_site': 'ProteinDomain',
            'Binding_site': 'ProteinDomain',
            'Conserved_site': 'ProteinDomain',
            'PTM': 'ProteinDomain',
        }

        for ipr_id, data in self.entries.items():
            entry_type = data['type']
            # Use a unified label for all InterPro entries
            label = "InterProEntry"

            props = {
                'name': data['name'],
                'short_name': data['short_name'],
                'entry_type': entry_type,
                'go_terms': data['go_terms'],
            }

            yield (ipr_id, label, props)
            node_count += 1

        logger.info(f"InterPro: Generated {node_count} InterProEntry nodes")

    def get_edges(self):
        """
        Generate hierarchy edges between InterPro entries.
        Yields: (id, source, target, label, properties)
        """
        logger.info("InterPro: Generating edges...")
        hier_count = 0

        for parent_id, child_id in self.hierarchy:
            yield (
                None,
                child_id,
                parent_id,
                "InterProChildOf",
                {}
            )
            hier_count += 1

        logger.info(f"InterPro: Generated {hier_count} InterProChildOf edges")
