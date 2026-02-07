"""
Rhea (Biochemical Reactions) Adapter for BioCypher.

Loads reaction data from Rhea TSV files and generates:
- BiochemicalReaction nodes (master + directional variants)
- ReactionVariant edges (variants → master)
- HasSubstrate/HasProduct edges (reactions → ChEBI compounds) via RDF parsing
"""

import gzip
import json
import pandas as pd
from biocypher._logger import logger


class RheaAdapter:
    def __init__(self, data_dir="template_package/data/rhea"):
        self.data_dir = data_dir
        self.directions = None
        self.ec_map = {}
        self.xrefs_map = {}
        self.smiles_map = {}
        self.chebi_names = {}
        self.participants = {}  # reaction_id -> {substrates: [], products: []}
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
        """Load all Rhea TSV files."""
        # 1. Load directions (master -> LR, RL, BI)
        logger.info("Rhea: Loading directions...")
        self.directions = pd.read_csv(
            f"{self.data_dir}/rhea-directions.tsv",
            sep='\t', dtype=str
        )
        logger.info(f"Rhea: Loaded {len(self.directions)} master reactions")

        # 2. Load EC number mappings
        logger.info("Rhea: Loading EC mappings...")
        ec_df = pd.read_csv(
            f"{self.data_dir}/rhea2ec.tsv",
            sep='\t', dtype=str
        )
        for _, row in ec_df.iterrows():
            rid = str(row.get('RHEA_ID', ''))
            ec = str(row.get('ID', ''))
            if rid and ec:
                if rid not in self.ec_map:
                    self.ec_map[rid] = []
                self.ec_map[rid].append(ec)
        logger.info(f"Rhea: Loaded EC mappings for {len(self.ec_map)} reactions")

        # 3. Load cross-references
        logger.info("Rhea: Loading cross-references...")
        xref_df = pd.read_csv(
            f"{self.data_dir}/rhea2xrefs.tsv",
            sep='\t', dtype=str
        )
        for _, row in xref_df.iterrows():
            rid = str(row.get('RHEA_ID', ''))
            db = str(row.get('DB', ''))
            xid = str(row.get('ID', ''))
            if rid and db and xid:
                if rid not in self.xrefs_map:
                    self.xrefs_map[rid] = {}
                if db not in self.xrefs_map[rid]:
                    self.xrefs_map[rid][db] = []
                self.xrefs_map[rid][db].append(xid)
        logger.info(f"Rhea: Loaded xrefs for {len(self.xrefs_map)} reactions")

        # 4. Load SMILES
        logger.info("Rhea: Loading SMILES...")
        try:
            with open(f"{self.data_dir}/rhea-reaction-smiles.tsv", 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        self.smiles_map[parts[0]] = parts[1]
        except FileNotFoundError:
            logger.warning("Rhea: SMILES file not found, skipping")
        logger.info(f"Rhea: Loaded {len(self.smiles_map)} SMILES")

        # 5. Load ChEBI participant names
        logger.info("Rhea: Loading ChEBI names...")
        try:
            chebi_df = pd.read_csv(
                f"{self.data_dir}/chebiId_name.tsv",
                sep='\t', dtype=str, header=None,
                names=['chebi_id', 'name']
            )
            for _, row in chebi_df.iterrows():
                cid = str(row.get('chebi_id', ''))
                name = str(row.get('name', ''))
                if cid:
                    self.chebi_names[cid] = name
        except Exception as e:
            logger.warning(f"Rhea: Could not load ChEBI names: {e}")
        logger.info(f"Rhea: Loaded {len(self.chebi_names)} ChEBI participant names")

        # 6. Parse RDF for reaction participants (substrates/products)
        self._parse_rdf_participants()

    def _parse_rdf_participants(self):
        """Parse rhea.rdf.gz to extract reaction substrates and products."""
        logger.info("Rhea: Parsing RDF for reaction participants...")
        rdf_path = f"{self.data_dir}/rhea.rdf.gz"

        try:
            import gzip
            import rdflib
            g = rdflib.Graph()
            logger.info("Rhea: Loading RDF graph (decompressing gz first)...")
            with gzip.open(rdf_path, 'rb') as gz_file:
                g.parse(gz_file, format='xml')
            logger.info(f"Rhea: RDF graph loaded with {len(g)} triples")

            # Namespaces
            RH = rdflib.Namespace("http://rdf.rhea-db.org/")

            # Build compound_id -> ChEBI accession map
            # Compounds have rh:accession like "CHEBI:15377"
            compound_chebi = {}
            for compound, _, acc in g.triples((None, RH.accession, None)):
                compound_uri = str(compound)
                acc_str = str(acc)
                if 'Compound_' in compound_uri and acc_str.startswith('CHEBI:'):
                    compound_chebi[compound_uri] = acc_str

            logger.info(f"Rhea: Mapped {len(compound_chebi)} RDF compounds to ChEBI IDs")

            # Find all reactions and their sides
            reaction_sides = {}

            # Helper to get ChEBI IDs from a reaction side
            def get_side_compounds(side_uri):
                chebi_ids = []
                for _, _, participant in g.triples((side_uri, RH.contains, None)):
                    # Participant has rh:compound -> Compound_XXXX
                    for _, _, compound_ref in g.triples((participant, RH.compound, None)):
                        compound_uri = str(compound_ref)
                        chebi_id = compound_chebi.get(compound_uri)
                        if chebi_id:
                            chebi_ids.append(chebi_id)
                return chebi_ids

            # Get substrates (left side) for each reaction
            for rxn, _, side in g.triples((None, RH.substrates, None)):
                rxn_id = str(rxn).split('/')[-1]
                if rxn_id not in reaction_sides:
                    reaction_sides[rxn_id] = {'left': [], 'right': []}
                reaction_sides[rxn_id]['left'].extend(get_side_compounds(side))

            # Get products (right side)
            for rxn, _, side in g.triples((None, RH.products, None)):
                rxn_id = str(rxn).split('/')[-1]
                if rxn_id not in reaction_sides:
                    reaction_sides[rxn_id] = {'left': [], 'right': []}
                reaction_sides[rxn_id]['right'].extend(get_side_compounds(side))

            self.participants = reaction_sides
            logger.info(f"Rhea: Extracted participants for {len(reaction_sides)} reactions")

        except ImportError:
            logger.warning("Rhea: rdflib not available, skipping RDF parsing. "
                         "Install with: pip install rdflib")
        except Exception as e:
            logger.warning(f"Rhea: RDF parsing failed: {e}. "
                         "Continuing without substrate/product edges.")

    def get_nodes(self):
        """
        Generate BiochemicalReaction nodes.
        Creates 4 nodes per master reaction (UN, LR, RL, BI).
        Yields: (id, label, properties)
        """
        logger.info("Rhea: Generating reaction nodes...")
        node_count = 0

        for _, row in self.directions.iterrows():
            master_id = str(row.get('RHEA_ID_MASTER', ''))
            lr_id = str(row.get('RHEA_ID_LR', ''))
            rl_id = str(row.get('RHEA_ID_RL', ''))
            bi_id = str(row.get('RHEA_ID_BI', ''))

            variants = [
                (master_id, 'UN'),
                (lr_id, 'LR'),
                (rl_id, 'RL'),
                (bi_id, 'BI'),
            ]

            for rhea_id, direction in variants:
                if not rhea_id:
                    continue

                # EC numbers (usually on master/UN)
                ec_list = self.ec_map.get(rhea_id, [])

                # Cross-references
                xrefs = self.xrefs_map.get(rhea_id, {})

                # SMILES (only LR and RL)
                smiles = self._sanitize(self.smiles_map.get(rhea_id, ''))

                # Build equation from participants
                parts = self.participants.get(rhea_id, {})
                left_names = [self.chebi_names.get(c, c) for c in parts.get('left', [])]
                right_names = [self.chebi_names.get(c, c) for c in parts.get('right', [])]
                equation = ""
                if left_names and right_names:
                    equation = self._sanitize(
                        " + ".join(left_names) + " = " + " + ".join(right_names)
                    )

                props = {
                    'rhea_id': f"RHEA:{rhea_id}",
                    'direction': direction,
                    'master_id': f"RHEA:{master_id}",
                    'ec_number': ec_list,
                    'smiles': smiles,
                    'equation': equation,
                    'xrefs': json.dumps(xrefs) if xrefs else "",
                }

                yield (f"RHEA:{rhea_id}", "BiochemicalReaction", props)
                node_count += 1

        logger.info(f"Rhea: Generated {node_count} BiochemicalReaction nodes")

    def get_edges(self):
        """
        Generate edges:
        1. ReactionVariant: LR/RL/BI → master (UN)
        2. HasSubstrate: reaction → ChEBI compound
        3. HasProduct: reaction → ChEBI compound
        Yields: (id, source, target, label, properties)
        """
        logger.info("Rhea: Generating edges...")
        variant_count = 0
        substrate_count = 0
        product_count = 0

        # 1. Variant edges
        for _, row in self.directions.iterrows():
            master_id = str(row.get('RHEA_ID_MASTER', ''))
            lr_id = str(row.get('RHEA_ID_LR', ''))
            rl_id = str(row.get('RHEA_ID_RL', ''))
            bi_id = str(row.get('RHEA_ID_BI', ''))

            for variant_id, variant_type in [(lr_id, 'LR'), (rl_id, 'RL'), (bi_id, 'BI')]:
                if variant_id and master_id:
                    yield (
                        None,
                        f"RHEA:{variant_id}",
                        f"RHEA:{master_id}",
                        "ReactionVariant",
                        {'variant_type': variant_type}
                    )
                    variant_count += 1

        # 2. Substrate and product edges (from RDF participants)
        for rxn_id, sides in self.participants.items():
            for chebi_id in sides.get('left', []):
                if chebi_id.startswith('CHEBI:'):
                    yield (
                        None,
                        f"RHEA:{rxn_id}",
                        chebi_id,
                        "HasSubstrate",
                        {}
                    )
                    substrate_count += 1

            for chebi_id in sides.get('right', []):
                if chebi_id.startswith('CHEBI:'):
                    yield (
                        None,
                        f"RHEA:{rxn_id}",
                        chebi_id,
                        "HasProduct",
                        {}
                    )
                    product_count += 1

        logger.info(f"Rhea: Generated {variant_count} ReactionVariant, "
                     f"{substrate_count} HasSubstrate, {product_count} HasProduct edges")
