"""
ChEBI (Chemical Entities of Biological Interest) Adapter for BioCypher.

Loads chemical compound data from ChEBI TSV files and generates:
- ChemicalSubstance nodes (compounds with properties)
- SubclassOf edges (ontological hierarchy)
- ChemicalRelation edges (other chemical relationships)
"""

import gzip
import json
import pandas as pd
from biocypher._logger import logger


class ChEBIAdapter:
    def __init__(self, data_dir="template_package/data/chebi"):
        self.data_dir = data_dir
        self.compounds = None
        self.chemical_data = None
        self.names = None
        self.relations = None
        self.db_accessions = None
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        # Remove control chars
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load all ChEBI TSV files."""
        logger.info("ChEBI: Loading compounds...")
        self.compounds = pd.read_csv(
            f"{self.data_dir}/compounds.tsv.gz",
            sep='\t', compression='gzip',
            dtype=str, na_values=['\\N', '']
        )
        logger.info(f"ChEBI: Loaded {len(self.compounds)} compounds")

        logger.info("ChEBI: Loading chemical data...")
        self.chemical_data = pd.read_csv(
            f"{self.data_dir}/chemical_data.tsv.gz",
            sep='\t', compression='gzip',
            dtype=str, na_values=['\\N', '']
        )
        logger.info(f"ChEBI: Loaded {len(self.chemical_data)} chemical data records")

        logger.info("ChEBI: Loading names...")
        self.names = pd.read_csv(
            f"{self.data_dir}/names.tsv.gz",
            sep='\t', compression='gzip',
            dtype=str, na_values=['\\N', '']
        )
        logger.info(f"ChEBI: Loaded {len(self.names)} name records")

        logger.info("ChEBI: Loading database accessions...")
        self.db_accessions = pd.read_csv(
            f"{self.data_dir}/database_accession.tsv.gz",
            sep='\t', compression='gzip',
            dtype=str, na_values=['\\N', '']
        )
        logger.info(f"ChEBI: Loaded {len(self.db_accessions)} accession records")

        logger.info("ChEBI: Loading relations...")
        self.relations = pd.read_csv(
            f"{self.data_dir}/relation.tsv.gz",
            sep='\t', compression='gzip',
            dtype=str, na_values=['\\N', '']
        )
        logger.info(f"ChEBI: Loaded {len(self.relations)} relations")

        # Build lookup indices
        self._build_indices()

    def _build_indices(self):
        """Build lookup dictionaries for efficient property assembly."""
        # Map internal compound_id -> chebi_accession
        self.id_to_chebi = {}
        for _, row in self.compounds.iterrows():
            if pd.notna(row.get('id')) and pd.notna(row.get('chebi_accession')):
                self.id_to_chebi[str(row['id'])] = str(row['chebi_accession'])

        # Chemical data by compound_id
        self.chem_by_compound = {}
        for _, row in self.chemical_data.iterrows():
            cid = str(row.get('compound_id', ''))
            if cid:
                self.chem_by_compound[cid] = {
                    'formula': row.get('formula'),
                    'charge': row.get('charge'),
                    'mass': row.get('mass'),
                    'monoisotopic_mass': row.get('monoisotopic_mass'),
                }

        # Names by compound_id (collect synonyms)
        self.names_by_compound = {}
        for _, row in self.names.iterrows():
            cid = str(row.get('compound_id', ''))
            if cid and pd.notna(row.get('name')):
                if cid not in self.names_by_compound:
                    self.names_by_compound[cid] = []
                self.names_by_compound[cid].append(self._sanitize(str(row['name'])))

        # DB accessions by compound_id
        self.xrefs_by_compound = {}
        for _, row in self.db_accessions.iterrows():
            cid = str(row.get('compound_id', ''))
            if cid and pd.notna(row.get('accession_number')):
                if cid not in self.xrefs_by_compound:
                    self.xrefs_by_compound[cid] = {}
                db_type = str(row.get('type', 'UNKNOWN'))
                acc = str(row['accession_number'])
                if db_type not in self.xrefs_by_compound[cid]:
                    self.xrefs_by_compound[cid][db_type] = []
                self.xrefs_by_compound[cid][db_type].append(acc)

        logger.info(f"ChEBI: Built indices - {len(self.id_to_chebi)} compounds, "
                     f"{len(self.chem_by_compound)} chemical data, "
                     f"{len(self.names_by_compound)} name sets, "
                     f"{len(self.xrefs_by_compound)} xref sets")

    def get_nodes(self):
        """
        Generate ChemicalSubstance nodes from ChEBI compounds.
        Yields: (id, label, properties)
        """
        logger.info("ChEBI: Generating nodes...")
        node_count = 0
        skipped = 0

        for _, row in self.compounds.iterrows():
            # Get ChEBI accession as the node ID
            chebi_acc = row.get('chebi_accession')
            if pd.isna(chebi_acc):
                skipped += 1
                continue

            chebi_acc = str(chebi_acc)
            internal_id = str(row.get('id', ''))
            status = str(row.get('status_id', ''))

            # Skip obsolete/merged entries (status_id=1 is active/checked)
            # status_id: 1=active, 2=checked, 3=submitted, etc.
            # Keep all entries that have a chebi_accession

            # Build properties
            name = self._sanitize(row.get('name', ''))
            definition = self._sanitize(row.get('definition', ''))
            stars = row.get('stars')

            # Add chemical data if available
            chem = self.chem_by_compound.get(internal_id, {})
            formula = self._sanitize(chem.get('formula', ''))
            charge_str = chem.get('charge')
            mass_str = chem.get('mass')

            charge = None
            if charge_str and pd.notna(charge_str):
                try:
                    charge = int(float(charge_str))
                except (ValueError, TypeError):
                    charge = None

            mass = None
            if mass_str and pd.notna(mass_str):
                try:
                    mass = float(mass_str)
                except (ValueError, TypeError):
                    mass = None

            stars_int = None
            if stars and pd.notna(stars):
                try:
                    stars_int = int(float(stars))
                except (ValueError, TypeError):
                    stars_int = None

            # Collect synonyms
            synonyms = self.names_by_compound.get(internal_id, [])

            # Collect cross-references as JSON string
            xrefs = self.xrefs_by_compound.get(internal_id, {})
            xrefs_json = json.dumps(xrefs) if xrefs else ""

            props = {
                'name': name,
                'definition': definition,
                'formula': formula,
                'charge': charge,
                'mass': mass,
                'stars': stars_int,
                'synonyms': synonyms,
                'xrefs': xrefs_json,
            }

            yield (chebi_acc, "ChemicalSubstance", props)
            node_count += 1

        logger.info(f"ChEBI: Generated {node_count} ChemicalSubstance nodes (skipped {skipped})")

    def get_edges(self):
        """
        Generate edges from ChEBI relations.
        Yields: (id, source, target, label, properties)
        """
        logger.info("ChEBI: Generating edges...")
        subclass_count = 0
        relation_count = 0
        skipped = 0

        for _, row in self.relations.iterrows():
            init_id = str(row.get('init_id', ''))
            final_id = str(row.get('final_id', ''))
            rel_type = str(row.get('relation_type_id', ''))

            # Map internal IDs to ChEBI accessions
            source_chebi = self.id_to_chebi.get(init_id)
            target_chebi = self.id_to_chebi.get(final_id)

            if not source_chebi or not target_chebi:
                skipped += 1
                continue

            # Type 5 = "is_a" (subclass)
            if rel_type == '5':
                yield (
                    None,
                    source_chebi,
                    target_chebi,
                    "SubclassOf",
                    {}
                )
                subclass_count += 1
            else:
                # Map relation type IDs to names
                rel_names = {
                    '1': 'is_conjugate_acid_of',
                    '2': 'is_conjugate_base_of',
                    '3': 'is_tautomer_of',
                    '4': 'is_enantiomer_of',
                    '6': 'has_part',
                    '7': 'has_role',
                    '8': 'has_parent_hydride',
                    '9': 'is_substituent_group_from',
                    '10': 'has_functional_parent',
                }
                rel_name = rel_names.get(rel_type, f'relation_{rel_type}')

                yield (
                    None,
                    source_chebi,
                    target_chebi,
                    "ChemicalRelation",
                    {'relation_type': rel_name}
                )
                relation_count += 1

        logger.info(f"ChEBI: Generated {subclass_count} SubclassOf edges, "
                     f"{relation_count} ChemicalRelation edges (skipped {skipped})")
