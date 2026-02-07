"""
GeneNetwork Adapter for BioCypher.

Loads expression trait and dataset metadata from GeneNetwork
(genenetwork.org) and generates:
- ExpressionQTL edges linking gene probesets to their peak QTL loci

GeneNetwork provides systems genetics data including expression QTL
(eQTL) mapping results.  The traits JSON files contain per-probeset
records with LRS scores, peak loci, chromosomal positions, and gene
annotations.

Expected data directory: template_package/data/genenetwork/
Expected files: traits_*.json, datasets_*.json, api_species.json,
                sample_data_*.json
"""

import json
from pathlib import Path
from biocypher._logger import logger


class GeneNetworkAdapter:
    def __init__(self, data_dir="template_package/data/genenetwork"):
        self.data_dir = Path(data_dir)
        self.traits = []       # list of dicts from traits JSON
        self.datasets = []     # list of dicts from datasets JSON
        self.species = []      # list of dicts from species JSON
        self._load_data()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
        return text.strip()

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def _load_data(self):
        """Load GeneNetwork JSON data files."""
        if not self.data_dir.exists():
            logger.warning(
                f"GeneNetwork: Data directory not found: {self.data_dir}"
            )
            return

        # Load species
        species_file = self.data_dir / "api_species.json"
        if species_file.exists():
            try:
                with open(species_file, "r", encoding="utf-8") as f:
                    self.species = json.load(f)
                logger.info(
                    f"GeneNetwork: Loaded {len(self.species)} species"
                )
            except Exception as e:
                logger.warning(
                    f"GeneNetwork: Error reading species file: {e}"
                )

        # Load datasets
        for fpath in sorted(self.data_dir.glob("datasets_*.json")):
            try:
                with open(fpath, "r", encoding="utf-8") as f:
                    data = json.load(f)
                if isinstance(data, list):
                    self.datasets.extend(data)
                logger.info(
                    f"GeneNetwork: Loaded {len(data)} datasets "
                    f"from {fpath.name}"
                )
            except Exception as e:
                logger.warning(
                    f"GeneNetwork: Error reading {fpath.name}: {e}"
                )

        # Load traits (primary data)
        for fpath in sorted(self.data_dir.glob("traits_*.json")):
            try:
                with open(fpath, "r", encoding="utf-8") as f:
                    data = json.load(f)
                if isinstance(data, list):
                    self.traits.extend(data)
                logger.info(
                    f"GeneNetwork: Loaded {len(data)} traits "
                    f"from {fpath.name}"
                )
            except Exception as e:
                logger.warning(
                    f"GeneNetwork: Error reading {fpath.name}: {e}"
                )

        logger.info(
            f"GeneNetwork: Total loaded: {len(self.traits)} traits, "
            f"{len(self.datasets)} datasets"
        )

    # ------------------------------------------------------------------
    # Node generator
    # ------------------------------------------------------------------
    def get_nodes(self):
        """
        Generate gene expression probe nodes from GeneNetwork trait data.

        Each trait/probeset is represented as a node carrying its
        chromosomal location and mean expression value.

        Yields: (id, label, properties)
        """
        logger.info(
            f"GeneNetwork: Generating nodes from "
            f"{len(self.traits)} traits..."
        )
        count = 0

        for trait in self.traits:
            probe_name = trait.get("Name")
            if not probe_name:
                continue

            gene_symbol = trait.get("Symbol", "")
            description = trait.get("Description", "")
            chrom = trait.get("Chr", "")
            mb = trait.get("Mb")
            mean_expr = trait.get("Mean")

            props = {
                "gene_symbol": self._sanitize(gene_symbol),
                "description": self._sanitize(description),
                "chromosome": self._sanitize(str(chrom)),
                "source": "GeneNetwork",
            }
            if mb is not None:
                props["position_mb"] = mb
            if mean_expr is not None:
                props["mean_expression"] = mean_expr
            if trait.get("Aliases"):
                props["aliases"] = self._sanitize(trait["Aliases"])

            yield (
                f"genenetwork:{self._sanitize(probe_name)}",
                "gene expression probe",
                props,
            )
            count += 1

        logger.info(f"GeneNetwork: Generated {count} probe nodes")

    # ------------------------------------------------------------------
    # Edge generator
    # ------------------------------------------------------------------
    def get_edges(self):
        """
        Generate eQTL edges linking probesets to their peak QTL loci.

        Only traits with a non-empty Locus field are included.

        Yields: (id, source, target, label, properties)
        """
        logger.info(
            f"GeneNetwork: Generating eQTL edges from "
            f"{len(self.traits)} traits..."
        )
        count = 0

        for trait in self.traits:
            probe_name = trait.get("Name")
            locus = trait.get("Locus")
            if not probe_name or not locus:
                continue

            lrs = trait.get("LRS")
            additive = trait.get("Additive")
            p_value = trait.get("P-Value")
            peak_chr = trait.get("Peak Chr")
            peak_mb = trait.get("Peak Mb")

            props = {
                "source": "GeneNetwork",
            }
            if lrs is not None:
                props["lrs_score"] = lrs
            if additive is not None:
                props["additive_effect"] = additive
            if p_value is not None:
                props["p_value"] = p_value
            if peak_chr is not None:
                props["peak_chromosome"] = self._sanitize(str(peak_chr))
            if peak_mb is not None:
                props["peak_position_mb"] = peak_mb

            yield (
                None,
                f"genenetwork:{self._sanitize(probe_name)}",
                f"genenetwork_locus:{self._sanitize(locus)}",
                "expression quantitative trait locus",
                props,
            )
            count += 1

        logger.info(f"GeneNetwork: Generated {count} eQTL edges")
