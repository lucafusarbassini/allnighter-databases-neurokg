"""
Brain RNA-Seq Adapter for BioCypher.

Loads cell-type-specific gene expression data from the Brain RNA-Seq
project (brainrnaseq.org) and generates:
- CellTypeExpression edges linking genes to cell types

The data comprises two CSV datasets:
- fe-wp-dataset-120.csv  (mouse, Mus musculus)
- fe-wp-dataset-124.csv  (human, Homo sapiens)

Each CSV has a gene_id column, an id column, and multiple replicate
columns per cell type.  Expression values are averaged across replicates
for each cell type.

Expected data directory: template_package/data/brain_rnaseq/
Expected files: fe-wp-dataset-120.csv, fe-wp-dataset-124.csv
"""

import csv
import re
from pathlib import Path
from biocypher._logger import logger


# Map filename to species
_FILE_SPECIES = {
    "fe-wp-dataset-120.csv": "Mus musculus",
    "fe-wp-dataset-124.csv": "Homo sapiens",
}


class BrainRNASeqAdapter:
    def __init__(self, data_dir="template_package/data/brain_rnaseq"):
        self.data_dir = Path(data_dir)
        self.expression_records = []  # list of dicts
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

    @staticmethod
    def _extract_cell_types(headers):
        """
        Group replicate columns by cell type.

        Column names follow the pattern ``<cell_type>_<replicate_number>``.
        Returns a dict mapping cell_type_name -> [column_name, ...].
        """
        cell_types = {}
        skip = {"gene_id", "id"}
        for h in headers:
            if h in skip:
                continue
            # Strip trailing _<digits> to get the base cell type name
            # Also handle columns like "oligodendrocytes_average_count"
            # that are not replicates -- include them as single-element
            m = re.match(r"^(.+?)_(\d+)$", h)
            if m:
                base = m.group(1)
            else:
                base = h
            cell_types.setdefault(base, []).append(h)
        return cell_types

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def _load_data(self):
        """Load CSV expression data files."""
        if not self.data_dir.exists():
            logger.warning(
                f"BrainRNASeq: Data directory not found: {self.data_dir}"
            )
            return

        csv_files = sorted(self.data_dir.glob("*.csv"))
        if not csv_files:
            logger.warning(
                f"BrainRNASeq: No CSV files found in {self.data_dir}"
            )
            return

        for fpath in csv_files:
            species = _FILE_SPECIES.get(fpath.name, "unknown")
            try:
                with open(fpath, "r", encoding="utf-8", errors="replace") as f:
                    reader = csv.DictReader(f)
                    headers = reader.fieldnames
                    if not headers:
                        logger.warning(
                            f"BrainRNASeq: Empty CSV {fpath.name}, skipping"
                        )
                        continue

                    cell_type_cols = self._extract_cell_types(headers)

                    for row in reader:
                        raw_gene_id = row.get("gene_id", "")
                        row_id = row.get("id", "")

                        # Extract gene symbol from strings like
                        # "A1BG - Homo sapiens" or "0610005C13Rik - Mus musculus"
                        gene_symbol = raw_gene_id.split(" - ")[0].strip()
                        if not gene_symbol:
                            continue

                        for cell_type, columns in cell_type_cols.items():
                            # Skip non-numeric aggregate columns
                            if cell_type in (
                                "oligodendrocytes_average_count",
                                "oligodendrocytes_standard_deviation",
                            ):
                                continue

                            values = []
                            for col in columns:
                                try:
                                    v = float(row.get(col, ""))
                                    values.append(v)
                                except (ValueError, TypeError):
                                    pass

                            if not values:
                                continue

                            avg_expr = sum(values) / len(values)

                            self.expression_records.append({
                                "gene_symbol": gene_symbol,
                                "gene_raw_id": row_id,
                                "cell_type": cell_type,
                                "expression_level": round(avg_expr, 6),
                                "n_replicates": len(values),
                                "species": species,
                            })

            except Exception as e:
                logger.warning(
                    f"BrainRNASeq: Error reading {fpath.name}: {e}"
                )

        logger.info(
            f"BrainRNASeq: Loaded {len(self.expression_records)} "
            f"gene-cell-type expression records"
        )

    # ------------------------------------------------------------------
    # Node generator
    # ------------------------------------------------------------------
    def get_nodes(self):
        """
        Brain RNA-Seq does not introduce new node types.  Gene and cell-type
        nodes are expected to exist from other adapters (e.g. HGNC, Cell
        Ontology).  This method yields nothing.
        """
        logger.info(
            "BrainRNASeq: No new nodes (uses existing Gene / CellType nodes)"
        )
        return
        yield  # makes this a generator

    # ------------------------------------------------------------------
    # Edge generator
    # ------------------------------------------------------------------
    def get_edges(self):
        """
        Generate CellTypeExpression edges (gene expressed in cell type).

        Yields: (id, source, target, label, properties)
        """
        logger.info(
            f"BrainRNASeq: Generating edges from "
            f"{len(self.expression_records)} expression records..."
        )
        count = 0

        for rec in self.expression_records:
            gene = self._sanitize(rec["gene_symbol"])
            cell_type = self._sanitize(
                rec["cell_type"].replace("_", " ")
            )

            props = {
                "expression_level": rec["expression_level"],
                "species": [rec["species"]],
                "cell_type": cell_type,
                "n_replicates": rec["n_replicates"],
                "source": "BrainRNASeq",
            }

            yield (
                None,
                gene,
                f"brain_rnaseq:{rec['cell_type']}",
                "gene expressed in cell type",
                props,
            )
            count += 1

        logger.info(
            f"BrainRNASeq: Generated {count} CellTypeExpression edges"
        )
