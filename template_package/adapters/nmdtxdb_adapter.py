"""
NMDtxDB Adapter for BioCypher.

Loads transcript-level NMD (nonsense-mediated mRNA decay) annotations from
the NMDtxDB database and generates:
- Transcript nodes with NMD annotation properties
- GeneHasTranscript edges linking genes to their transcripts

The upstream data file is an R serialized object (database.RDS).  If the
optional ``pyreadr`` package is available the adapter will parse it
natively; otherwise it emits a warning and produces no records.

Expected data directory: template_package/data/nmdtxdb/
Expected file: database.RDS
"""

from pathlib import Path
from biocypher._logger import logger


class NMDtxDBAdapter:
    def __init__(self, data_dir="template_package/data/nmdtxdb"):
        self.data_dir = Path(data_dir)
        self.transcripts = []  # list of dicts
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
        """Load NMDtxDB data from RDS file using pyreadr if available."""
        rds_path = self.data_dir / "database.RDS"
        if not rds_path.exists():
            logger.warning(
                f"NMDtxDB: RDS file not found at {rds_path}"
            )
            return

        try:
            import pyreadr
        except ImportError:
            logger.warning(
                "NMDtxDB: 'pyreadr' package is not installed. "
                "Cannot parse database.RDS. Install with: "
                "pip install pyreadr"
            )
            return

        try:
            result = pyreadr.read_r(str(rds_path))
            # read_r returns an OrderedDict of DataFrames; RDS has one key
            if not result:
                logger.warning("NMDtxDB: RDS file was empty")
                return

            # Get the first (and typically only) DataFrame
            df_key = list(result.keys())[0]
            df = result[df_key]

            logger.info(
                f"NMDtxDB: Read dataframe '{df_key}' with "
                f"{len(df)} rows, columns: {list(df.columns)}"
            )

            # Build transcript records from dataframe rows
            for _, row in df.iterrows():
                record = {}
                for col in df.columns:
                    record[col.lower().replace(" ", "_")] = row[col]
                self.transcripts.append(record)

        except Exception as e:
            logger.warning(f"NMDtxDB: Error reading RDS file: {e}")
            return

        logger.info(
            f"NMDtxDB: Loaded {len(self.transcripts)} transcript records"
        )

    # ------------------------------------------------------------------
    # Node generator
    # ------------------------------------------------------------------
    def get_nodes(self):
        """
        Generate Transcript nodes from NMDtxDB data.

        Each node carries NMD annotation properties such as PTC status,
        NMD sensitivity, and transcript biotype when available.

        Yields: (id, label, properties)
        """
        logger.info(
            f"NMDtxDB: Generating nodes from "
            f"{len(self.transcripts)} transcript records..."
        )
        count = 0

        for rec in self.transcripts:
            # Try common column names for transcript ID
            tx_id = (
                rec.get("transcript_id")
                or rec.get("tx_id")
                or rec.get("ensembl_transcript_id")
                or rec.get("name")
            )
            if not tx_id:
                continue

            tx_id = self._sanitize(str(tx_id))

            props = {
                "source": "NMDtxDB",
            }

            # Map known annotation fields
            field_map = {
                "gene_id": "gene_id",
                "gene_name": "gene_name",
                "gene_symbol": "gene_symbol",
                "tx_biotype": "transcript_biotype",
                "transcript_biotype": "transcript_biotype",
                "ptc_status": "ptc_status",
                "ptc": "ptc_status",
                "nmd": "nmd_sensitive",
                "nmd_sensitive": "nmd_sensitive",
                "nmd_status": "nmd_status",
                "chr": "chromosome",
                "chromosome": "chromosome",
                "start": "start",
                "end": "end",
                "strand": "strand",
            }

            for src_key, dst_key in field_map.items():
                val = rec.get(src_key)
                if val is not None:
                    props[dst_key] = self._sanitize(str(val))

            yield (tx_id, "transcript", props)
            count += 1

        logger.info(f"NMDtxDB: Generated {count} Transcript nodes")

    # ------------------------------------------------------------------
    # Edge generator
    # ------------------------------------------------------------------
    def get_edges(self):
        """
        Generate GeneHasTranscript edges linking gene IDs to transcripts.

        Yields: (id, source, target, label, properties)
        """
        logger.info(
            f"NMDtxDB: Generating edges from "
            f"{len(self.transcripts)} transcript records..."
        )
        count = 0

        for rec in self.transcripts:
            tx_id = (
                rec.get("transcript_id")
                or rec.get("tx_id")
                or rec.get("ensembl_transcript_id")
                or rec.get("name")
            )
            gene_id = (
                rec.get("gene_id")
                or rec.get("ensembl_gene_id")
                or rec.get("gene_name")
                or rec.get("gene_symbol")
            )

            if not tx_id or not gene_id:
                continue

            tx_id = self._sanitize(str(tx_id))
            gene_id = self._sanitize(str(gene_id))

            props = {
                "source": "NMDtxDB",
            }

            # Carry forward NMD annotation on the edge as well
            for key in ("ptc_status", "ptc", "nmd", "nmd_sensitive",
                        "nmd_status"):
                val = rec.get(key)
                if val is not None:
                    props[key] = self._sanitize(str(val))

            yield (None, gene_id, tx_id, "gene has transcript", props)
            count += 1

        logger.info(f"NMDtxDB: Generated {count} GeneHasTranscript edges")
