import pandas as pd
import biocypher
from biocypher._logger import logger

class LianaAdapter:
    def __init__(self):
        self.human_file = "template_package/data/liana_humanconsensus_db.parquet"
        self.mouse_file = "template_package/data/liana_mouseconsensus_db.parquet" # this already uses human id (ortholog) for mouse genes
        self.data = self._load_data()

    def _load_data(self):
        """
        Loads and merges the mouse and human parquet files.
        The parquet files are expected to have columns: 'source', 'target', 'species'
        """
        logger.info("Loading Parquet files...")

        try:
            df_human = pd.read_parquet(self.human_file)
            if "species" not in df_human.columns:
                df_human["species"] = "Homo sapiens"
        except FileNotFoundError:
            logger.warning(f"File {self.human_file} not found. Skipping.")
            df_human = pd.DataFrame()

        try:
            df_mouse = pd.read_parquet(self.mouse_file)
            if "species" not in df_mouse.columns:
                df_mouse["species"] = "Mus musculus"
        except FileNotFoundError:
            logger.warning(f"File {self.mouse_file} not found. Skipping.")
            df_mouse = pd.DataFrame()

        # Combine datasets
        full_df = pd.concat([df_human, df_mouse], ignore_index=True)
        merged_species = (
            full_df.groupby(["source", "target"])["species"]
              .apply(lambda x: sorted(set(x)))   # deduplicate + sort
              .reset_index()
        )

        # Validate required columns exist
        required_cols = ['source', 'target', 'species']
        missing_cols = [col for col in required_cols if col not in merged_species.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}. Found columns: {merged_species.columns.tolist()}")

        logger.info(f"Loaded {len(merged_species)} interactions from {full_df['species'].nunique()} species.")
        logger.info(f"Species distribution: only human -> {(merged_species['species'].apply(lambda x: x == ['Homo sapiens'])).sum()}, only mouse -> {(merged_species['species'].apply(lambda x: x == ['Mus musculus'])).sum()}, both -> {(merged_species['species'].apply(lambda x: x == ['Homo sapiens', 'Mus musculus'])).sum()}")
        
        return merged_species

    def get_nodes(self):
        """
        Generator yielding nodes.
        Format: (id, label, properties)
        
        Each unique protein (from ligand or receptor columns) becomes a node.
        If an entry starts with "COMPLEX:" (e.g. "COMPLEX:P04626_P21860") it is treated as a QuaternaryStructure.
        No species in nodes -> species is stored as a property in edges.
        """
        logger.info("Generating Nodes...")

        # Map id -> {"species": species, "label": "Gene"|"QuaternaryStructure"}
        node_to_label = {}

        for _, row in self.data.iterrows():
            for col in ("source", "target"):
                pid = row[col]

                if pd.isna(pid):
                    continue

                pid_str = str(pid)

                if pid_str.startswith("COMPLEX:") or "_" in pid_str:
                    label = "QuaternaryStructure"
                    pid_str = pid_str.replace("COMPLEX:", "")
                else:
                    label = "Gene"

                if pid_str not in node_to_label:
                    node_to_label[pid_str] = label

        # Yield nodes
        for protein_id, label in node_to_label.items():
            yield (
                protein_id,
                label,
                {}
            )

        logger.info(f"Generated {len(node_to_label)} unique protein/complex nodes.")

    def get_edges(self):
        """
        Generator yielding edges.
        Format: (id, source_id, target_id, label, properties)
        
        Each row represents a ligand-receptor interaction.
        """
        logger.info("Generating Edges...")
        for _, row in self.data.iterrows():
            ligand = row["source"]
            receptor = row["target"]
            species = row["species"]

            if ligand.startswith("COMPLEX:") or "_" in ligand:
                ligand = ligand.replace("COMPLEX:", "")
            if receptor.startswith("COMPLEX:") or "_" in receptor:
                receptor = receptor.replace("COMPLEX:", "")

            # Additional properties beyond species can be added here
            properties = {"species": species}
            #rel_id = f'{ligand}_{receptor}_{len(species)}'  # Unique edge id
            # With id = None, Biocypher dedupes based on ONLY source, target, label ! (not properties) -> keeps the first occurrence of the edge.  
            # Okay here since we merged same source-target pair with different species into one with a list of species.
            yield (
                None,  
                ligand,      # Source (ligand)      
                receptor,    # Target (receptor)
                "LigandReceptorInteraction",
                properties
            )
        
        logger.info(f"Generated {len(self.data)} interaction edges.")

