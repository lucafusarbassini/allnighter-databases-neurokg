"""
LIPID MAPS Structure Database Adapter for BioCypher.

Loads lipid structure data from LIPIDMAPS SDF files and generates:
- Lipid nodes (chemical structures with classification)
- LipidCategory nodes (hierarchical classification)
- LipidClassifiedAs edges (lipid → category)
- EquivalentTo edges (lipid → ChEBI compound)
"""

import re
import json
from pathlib import Path
from biocypher._logger import logger


class LIPIDMAPSAdapter:
    def __init__(self, data_dir="template_package/data/lipidmaps"):
        self.data_dir = data_dir
        self.lipids = {}
        self.categories = {}
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _clean_value(self, val):
        """Clean a value: return empty string for missing/dash values."""
        if val is None or val == '-' or val == '':
            return ''
        return self._sanitize(val)

    def _load_data(self):
        """Load LIPIDMAPS SDF file."""
        sdf_path = f"{self.data_dir}/structures_extended.sdf"
        if not Path(sdf_path).exists():
            sdf_path = f"{self.data_dir}/structures.sdf"

        logger.info(f"LIPIDMAPS: Parsing SDF file: {sdf_path}")
        self.lipids = self._parse_sdf(sdf_path)
        logger.info(f"LIPIDMAPS: Loaded {len(self.lipids)} lipids")

        # Build category hierarchy
        self._build_categories()

    def _parse_sdf(self, filepath):
        """Parse SDF file and extract metadata fields."""
        records = {}
        current_record = {}
        current_field = None
        in_structure = True
        lm_id = None
        line_num = 0

        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.rstrip('\n')
                line_num += 1

                if line == '$$$$':
                    # End of record
                    if lm_id and current_record:
                        records[lm_id] = current_record
                    current_record = {}
                    current_field = None
                    in_structure = True
                    lm_id = None

                    # Log progress every 10000 records
                    if len(records) % 10000 == 0 and len(records) > 0:
                        logger.info(f"LIPIDMAPS: Parsed {len(records)} records...")
                    continue

                # First line of a record is the molecule name (ID)
                if in_structure and lm_id is None and line.strip():
                    lm_id = line.strip()

                # Check if we're past the structure block (field header)
                if line.startswith('> <'):
                    in_structure = False
                    # Extract field name from "> <FIELD_NAME>"
                    current_field = line[3:]
                    if current_field.endswith('>'):
                        current_field = current_field[:-1]
                    current_record[current_field] = None
                    continue

                if in_structure:
                    continue

                # Parse field values
                if current_field is not None and line:
                    if current_record[current_field] is None:
                        current_record[current_field] = line
                    else:
                        current_record[current_field] += '\n' + line

        # Capture the last record
        if lm_id and current_record:
            records[lm_id] = current_record

        return records

    def _extract_category_code(self, category_str):
        """Extract category code from string like 'Fatty Acyls [FA]'."""
        if not category_str:
            return None
        match = re.search(r'\[([^\]]+)\]', category_str)
        if match:
            return match.group(1)
        return category_str

    def _build_categories(self):
        """Build category hierarchy from lipid data."""
        logger.info("LIPIDMAPS: Building category hierarchy...")
        categories = {}

        for lm_id, data in self.lipids.items():
            # Process each classification level
            levels = [
                ('CATEGORY', 'category'),
                ('MAIN_CLASS', 'main_class'),
                ('SUB_CLASS', 'sub_class'),
                ('CLASS_LEVEL4', 'level4'),
            ]

            for field, level_name in levels:
                val = data.get(field, '')
                if val and val != '-':
                    code = self._extract_category_code(val)
                    if code and code not in categories:
                        categories[code] = {
                            'code': code,
                            'name': self._sanitize(val),
                            'level': level_name,
                        }

        self.categories = categories
        logger.info(f"LIPIDMAPS: Built {len(self.categories)} category nodes")

    def _build_lipid_properties(self, data):
        """Build node properties from SDF record data."""
        props = {
            'name': self._clean_value(data.get('COMMON_NAME')),
            'systematic_name': self._clean_value(data.get('SYSTEMATIC_NAME')),
            'abbreviation': self._clean_value(data.get('ABBREVIATION')),
            'formula': self._clean_value(data.get('FORMULA')),
            'inchi': self._clean_value(data.get('INCHI')),
            'inchi_key': self._clean_value(data.get('INCHI_KEY')),
            'smiles': self._clean_value(data.get('SMILES')),
            'category': self._clean_value(data.get('CATEGORY')),
            'main_class': self._clean_value(data.get('MAIN_CLASS')),
            'sub_class': self._clean_value(data.get('SUB_CLASS')),
            'class_level4': self._clean_value(data.get('CLASS_LEVEL4')),
            'source': 'LIPIDMAPS',
        }

        # Mass
        mass_str = data.get('EXACT_MASS')
        if mass_str and mass_str != '-':
            try:
                props['mass'] = float(mass_str)
            except (ValueError, TypeError):
                props['mass'] = None
        else:
            props['mass'] = None

        # Synonyms
        synonyms = []
        syn_val = data.get('SYNONYMS', '')
        if syn_val and syn_val != '-':
            synonyms = [self._sanitize(s.strip()) for s in syn_val.split(';') if s.strip() and s.strip() != '-']
        props['synonyms'] = synonyms

        # Cross-references as JSON
        xrefs = {}
        xref_fields = [
            ('PUBCHEM_CID', 'PUBCHEM'),
            ('CHEBI_ID', 'CHEBI'),
            ('KEGG_ID', 'KEGG'),
            ('HMDB_ID', 'HMDB'),
            ('LIPIDBANK_ID', 'LIPIDBANK'),
            ('SWISSLIPIDS_ID', 'SWISSLIPIDS'),
            ('CAYMAN_ID', 'CAYMAN'),
        ]
        for field, key in xref_fields:
            val = data.get(field, '')
            if val and val != '-':
                xrefs[key] = val
        props['xrefs'] = json.dumps(xrefs) if xrefs else ""

        return props

    def get_nodes(self):
        """
        Generate Lipid and LipidCategory nodes.
        Yields: (id, label, properties)
        """
        logger.info("LIPIDMAPS: Generating nodes...")
        lipid_count = 0
        cat_count = 0

        # Generate Lipid nodes
        for lm_id, data in self.lipids.items():
            props = self._build_lipid_properties(data)
            yield (lm_id, "Lipid", props)
            lipid_count += 1

        # Generate LipidCategory nodes
        for cat_code, cat_data in self.categories.items():
            yield (cat_code, "LipidCategory", cat_data)
            cat_count += 1

        logger.info(f"LIPIDMAPS: Generated {lipid_count} Lipid nodes, {cat_count} LipidCategory nodes")

    def get_edges(self):
        """
        Generate classification and equivalence edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("LIPIDMAPS: Generating edges...")
        classify_count = 0
        equiv_count = 0

        for lm_id, data in self.lipids.items():
            # Classification edges
            levels = [
                ('CATEGORY', 'category'),
                ('MAIN_CLASS', 'main_class'),
                ('SUB_CLASS', 'sub_class'),
                ('CLASS_LEVEL4', 'level4'),
            ]

            for field, level_name in levels:
                val = data.get(field, '')
                if val and val != '-':
                    code = self._extract_category_code(val)
                    if code:
                        yield (
                            None,
                            lm_id,
                            code,
                            "LipidClassifiedAs",
                            {
                                'classification_level': level_name,
                                'classification_type': 'primary',
                            }
                        )
                        classify_count += 1

            # ChEBI equivalence edge
            chebi_id = data.get('CHEBI_ID', '')
            if chebi_id and chebi_id != '-':
                # Normalize to CHEBI:XXXXX format
                if not chebi_id.startswith('CHEBI:'):
                    chebi_id = f"CHEBI:{chebi_id}"
                yield (
                    None,
                    lm_id,
                    chebi_id,
                    "EquivalentTo",
                    {'source': 'LIPIDMAPS'}
                )
                equiv_count += 1

        logger.info(f"LIPIDMAPS: Generated {classify_count} LipidClassifiedAs, "
                     f"{equiv_count} EquivalentTo edges")
