"""
KEGG (Kyoto Encyclopedia of Genes and Genomes) Adapter for BioCypher.

Loads KEGG data from downloaded TSV/text files and generates:
- KEGGCompound nodes (metabolites)
- KEGGReaction nodes (metabolic reactions)
- KEGGPathway nodes (metabolic pathways)
- CompoundParticipatesInReaction edges
- ReactionInPathway edges
- KEGGCrossRef edges (KEGG compound → ChEBI)
"""

import json
import re
from pathlib import Path
from biocypher._logger import logger


class KEGGAdapter:
    def __init__(self, data_dir="template_package/data/kegg"):
        self.data_dir = data_dir
        self.compounds = {}
        self.reactions = {}
        self.pathways = {}
        self.compound_to_chebi = {}
        self.compound_to_pubchem = {}
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
        """Load all KEGG data files."""
        self._load_lists()
        self._load_xrefs()
        self._load_entries()

    def _load_lists(self):
        """Load KEGG list files (compounds, reactions, pathways)."""
        # Compounds
        logger.info("KEGG: Loading compound list...")
        comp_path = f"{self.data_dir}/lists/compound_list.tsv"
        if Path(comp_path).exists():
            with open(comp_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        cid = parts[0].strip()
                        names_str = parts[1].strip()
                        # First name before semicolons is the primary name
                        names = [n.strip() for n in names_str.split(';')]
                        self.compounds[cid] = {
                            'name': self._sanitize(names[0]) if names else '',
                            'synonyms': [self._sanitize(n) for n in names[1:] if n.strip()],
                        }
        logger.info(f"KEGG: Loaded {len(self.compounds)} compounds")

        # Reactions
        logger.info("KEGG: Loading reaction list...")
        rxn_path = f"{self.data_dir}/lists/reaction_list.tsv"
        if Path(rxn_path).exists():
            with open(rxn_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        rid = parts[0].strip()
                        desc = parts[1].strip()
                        # Split on semicolons - first part is enzyme name, second is equation
                        desc_parts = desc.split(';', 1)
                        enzyme_name = self._sanitize(desc_parts[0].strip()) if desc_parts else ''
                        equation = self._sanitize(desc_parts[1].strip()) if len(desc_parts) > 1 else ''
                        self.reactions[rid] = {
                            'enzyme_name': enzyme_name,
                            'equation': equation,
                        }
        logger.info(f"KEGG: Loaded {len(self.reactions)} reactions")

        # Pathways
        logger.info("KEGG: Loading pathway list...")
        pw_path = f"{self.data_dir}/lists/pathway_list.tsv"
        if Path(pw_path).exists():
            with open(pw_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        pid = parts[0].strip()
                        name = self._sanitize(parts[1].strip())
                        self.pathways[pid] = {
                            'name': name,
                        }
        logger.info(f"KEGG: Loaded {len(self.pathways)} pathways")

    def _load_xrefs(self):
        """Load cross-reference files."""
        # Compound to ChEBI
        chebi_path = f"{self.data_dir}/xrefs/compound_to_chebi.tsv"
        if Path(chebi_path).exists():
            with open(chebi_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        # cpd:C00001 -> chebi:15377
                        kegg_id = parts[0].replace('cpd:', '')
                        chebi_id = parts[1].replace('chebi:', 'CHEBI:')
                        self.compound_to_chebi[kegg_id] = chebi_id
        logger.info(f"KEGG: Loaded {len(self.compound_to_chebi)} compound-ChEBI mappings")

        # Compound to PubChem
        pubchem_path = f"{self.data_dir}/xrefs/compound_to_pubchem.tsv"
        if Path(pubchem_path).exists():
            with open(pubchem_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        kegg_id = parts[0].replace('cpd:', '')
                        pubchem_id = parts[1].replace('pubchem:', '')
                        self.compound_to_pubchem[kegg_id] = pubchem_id
        logger.info(f"KEGG: Loaded {len(self.compound_to_pubchem)} compound-PubChem mappings")

    def _load_entries(self):
        """Load detailed KEGG entry data from batch files."""
        # Parse compound entries
        compound_dir = Path(f"{self.data_dir}/entries/compound")
        if compound_dir.exists():
            logger.info("KEGG: Loading compound entries...")
            for batch_file in sorted(compound_dir.glob("batch_*.txt")):
                self._parse_kegg_flat_file(batch_file, 'compound')

        # Parse reaction entries
        reaction_dir = Path(f"{self.data_dir}/entries/reaction")
        if reaction_dir.exists():
            logger.info("KEGG: Loading reaction entries...")
            for batch_file in sorted(reaction_dir.glob("batch_*.txt")):
                self._parse_kegg_flat_file(batch_file, 'reaction')

        # Parse pathway entries
        pathway_dir = Path(f"{self.data_dir}/entries/pathway")
        if pathway_dir.exists():
            logger.info("KEGG: Loading pathway entries...")
            for batch_file in sorted(pathway_dir.glob("batch_*.txt")):
                self._parse_kegg_flat_file(batch_file, 'pathway')

    def _parse_kegg_flat_file(self, filepath, entry_type):
        """Parse a KEGG flat file format and merge data into existing records."""
        current_id = None
        current_field = None
        current_data = {}

        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.rstrip('\n')

                # End of entry
                if line.startswith('///'):
                    if current_id and current_data:
                        self._merge_entry(current_id, current_data, entry_type)
                    current_id = None
                    current_field = None
                    current_data = {}
                    continue

                # Field header (12-char fixed width)
                if line and not line[0].isspace() and len(line) > 12:
                    field_name = line[:12].strip()
                    value = line[12:].strip()

                    if field_name == 'ENTRY':
                        # Extract ID from entry line
                        parts = value.split()
                        current_id = parts[0] if parts else None
                        current_field = 'ENTRY'
                    elif field_name:
                        current_field = field_name
                        if current_field not in current_data:
                            current_data[current_field] = value
                        else:
                            current_data[current_field] += ' ' + value
                elif line.startswith(' ' * 12) and current_field:
                    # Continuation line
                    value = line[12:].strip()
                    if current_field in current_data:
                        current_data[current_field] += ' ' + value
                    else:
                        current_data[current_field] = value
                elif line.startswith('  ') and current_field:
                    # Alternative continuation format
                    value = line.strip()
                    if current_field in current_data:
                        current_data[current_field] += ' ' + value
                    else:
                        current_data[current_field] = value

        # Handle last entry in file
        if current_id and current_data:
            self._merge_entry(current_id, current_data, entry_type)

    def _merge_entry(self, entry_id, data, entry_type):
        """Merge parsed entry data into existing records."""
        if entry_type == 'compound':
            if entry_id in self.compounds:
                # Add formula, mass, etc.
                if 'FORMULA' in data:
                    self.compounds[entry_id]['formula'] = self._sanitize(data['FORMULA'])
                if 'EXACT_MASS' in data:
                    try:
                        self.compounds[entry_id]['mass'] = float(data['EXACT_MASS'])
                    except (ValueError, TypeError):
                        pass
                if 'MOL_WEIGHT' in data:
                    try:
                        self.compounds[entry_id]['mol_weight'] = float(data['MOL_WEIGHT'])
                    except (ValueError, TypeError):
                        pass
                if 'REACTION' in data:
                    rxn_ids = data['REACTION'].split()
                    self.compounds[entry_id]['reactions'] = [r for r in rxn_ids if r.startswith('R')]
                if 'PATHWAY' in data:
                    # Parse pathway references
                    pathways = re.findall(r'(map\d+)', data.get('PATHWAY', ''))
                    self.compounds[entry_id]['pathways'] = pathways

        elif entry_type == 'reaction':
            if entry_id in self.reactions:
                if 'EQUATION' in data:
                    self.reactions[entry_id]['equation'] = self._sanitize(data['EQUATION'])
                if 'ENZYME' in data:
                    enzymes = data['ENZYME'].split()
                    self.reactions[entry_id]['ec_numbers'] = enzymes
                if 'PATHWAY' in data:
                    pathways = re.findall(r'(rn\d+|map\d+)', data.get('PATHWAY', ''))
                    self.reactions[entry_id]['pathways'] = pathways
                # Extract compound participants from equation
                if 'EQUATION' in data:
                    compounds = re.findall(r'(C\d{5})', data['EQUATION'])
                    self.reactions[entry_id]['compounds'] = compounds

        elif entry_type == 'pathway':
            if entry_id in self.pathways:
                if 'DESCRIPTION' in data:
                    self.pathways[entry_id]['description'] = self._sanitize(data['DESCRIPTION'])
                if 'CLASS' in data:
                    self.pathways[entry_id]['pathway_class'] = self._sanitize(data['CLASS'])

    def get_nodes(self):
        """
        Generate KEGG compound, reaction, and pathway nodes.
        Yields: (id, label, properties)
        """
        logger.info("KEGG: Generating nodes...")
        comp_count = 0
        rxn_count = 0
        pw_count = 0

        # Compound nodes
        for cid, data in self.compounds.items():
            xrefs = {}
            if cid in self.compound_to_chebi:
                xrefs['CHEBI'] = self.compound_to_chebi[cid]
            if cid in self.compound_to_pubchem:
                xrefs['PUBCHEM'] = self.compound_to_pubchem[cid]

            props = {
                'name': data.get('name', ''),
                'synonyms': data.get('synonyms', []),
                'formula': data.get('formula', ''),
                'mass': data.get('mass'),
                'xrefs': json.dumps(xrefs) if xrefs else '',
                'source': 'KEGG',
            }
            yield (f"KEGG:{cid}", "KEGGCompound", props)
            comp_count += 1

        # Reaction nodes
        for rid, data in self.reactions.items():
            props = {
                'enzyme_name': data.get('enzyme_name', ''),
                'equation': data.get('equation', ''),
                'ec_numbers': data.get('ec_numbers', []),
                'source': 'KEGG',
            }
            yield (f"KEGG:{rid}", "KEGGReaction", props)
            rxn_count += 1

        # Pathway nodes
        for pid, data in self.pathways.items():
            props = {
                'name': data.get('name', ''),
                'description': data.get('description', ''),
                'pathway_class': data.get('pathway_class', ''),
                'source': 'KEGG',
            }
            yield (f"KEGG:{pid}", "KEGGPathway", props)
            pw_count += 1

        logger.info(f"KEGG: Generated {comp_count} compound, {rxn_count} reaction, "
                     f"{pw_count} pathway nodes")

    def get_edges(self):
        """
        Generate KEGG relationship edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("KEGG: Generating edges...")
        comp_rxn_count = 0
        rxn_pw_count = 0
        xref_count = 0

        # Compound → Reaction participation edges (from reaction equations)
        for rid, data in self.reactions.items():
            for cid in data.get('compounds', []):
                if cid in self.compounds:
                    yield (
                        None,
                        f"KEGG:{cid}",
                        f"KEGG:{rid}",
                        "CompoundParticipatesInReaction",
                        {}
                    )
                    comp_rxn_count += 1

        # Reaction → Pathway edges
        for rid, data in self.reactions.items():
            for pid in data.get('pathways', []):
                if pid in self.pathways:
                    yield (
                        None,
                        f"KEGG:{rid}",
                        f"KEGG:{pid}",
                        "ReactionInPathway",
                        {}
                    )
                    rxn_pw_count += 1

        # KEGG compound → ChEBI cross-reference edges
        for cid, chebi_id in self.compound_to_chebi.items():
            if cid in self.compounds:
                yield (
                    None,
                    f"KEGG:{cid}",
                    chebi_id,
                    "EquivalentTo",
                    {'source': 'KEGG'}
                )
                xref_count += 1

        logger.info(f"KEGG: Generated {comp_rxn_count} CompoundParticipatesInReaction, "
                     f"{rxn_pw_count} ReactionInPathway, {xref_count} EquivalentTo edges")
