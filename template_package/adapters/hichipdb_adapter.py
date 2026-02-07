"""
HiChIPdb (3D Regulatory Interactions) Adapter for BioCypher.

Parses BEDPE chromatin loop files and ENCODE metadata to yield
ChromatinInteraction edges between two genomic loci.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class HiChIPdbAdapter:
    def __init__(self, data_dir="template_package/data/hichipdb"):
        self.data_dir = Path(data_dir)
        self.loops = []
        self.metadata = {}       # accession -> metadata dict
        self.fourdn_info = []    # 4DN experiment records
        self._load_data()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _sanitize(text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
        return text.strip()

    @staticmethod
    def _locus_id(chrom, start, end):
        """Canonical locus identifier."""
        return f"{chrom}:{start}-{end}"

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("HiChIPdb: data directory not found")
            return
        self._load_metadata()
        self._load_fourdn()
        self._load_bedpe_files()
        logger.info(
            f"HiChIPdb: Loaded {len(self.loops)} loops, "
            f"{len(self.metadata)} ENCODE metadata entries, "
            f"{len(self.fourdn_info)} 4DN experiments"
        )

    def _load_metadata(self):
        """Parse encode_loops_metadata.tsv."""
        fpath = self.data_dir / "encode_loops_metadata.tsv"
        if not fpath.exists():
            return
        with open(fpath, "r", errors="replace") as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if header is None:
                    header = parts
                    continue
                row = dict(zip(header, parts))
                acc = row.get("file_accession", "")
                if acc:
                    self.metadata[acc] = row

    def _load_fourdn(self):
        """Parse fourdn_hichip_experiments.json."""
        fpath = self.data_dir / "fourdn_hichip_experiments.json"
        if not fpath.exists():
            return
        with open(fpath, "r", errors="replace") as fh:
            data = json.load(fh)
            if isinstance(data, list):
                self.fourdn_info = data
            elif isinstance(data, dict):
                self.fourdn_info = data.get("results", data.get("@graph", [data]))

    def _load_bedpe_files(self):
        """Parse all *.bedpe files in the data directory."""
        for fpath in sorted(self.data_dir.glob("*.bedpe")):
            try:
                self._parse_bedpe(fpath)
            except Exception as exc:
                logger.warning(f"HiChIPdb: Error reading {fpath}: {exc}")

    def _parse_bedpe(self, fpath):
        """
        BEDPE columns:
          chr1  x1  x2  chr2  y1  y2  name  score  strand1  strand2
          color  observed  expectedBL  expectedDonut  expectedH  expectedV
          fdrBL  fdrDonut  fdrH  fdrV  numCollapsed
          centroid1  centroid2  radius
          highRes_start_1 highRes_end_1 highRes_start_2 highRes_end_2 ...
        """
        # Derive the ENCODE accession from filename (e.g. ENCFF661SAZ)
        accession = fpath.stem.replace("_loops", "")

        with open(fpath, "r", errors="replace") as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                # skip comment lines (juicer version etc.)
                if line.startswith("#") and header is not None:
                    continue
                parts = line.split("\t")
                if header is None:
                    # first line that starts with # is the header
                    if line.startswith("#"):
                        parts[0] = parts[0].lstrip("#")
                    header = parts
                    continue
                if len(parts) < 6:
                    continue
                row = dict(zip(header, parts))

                chr1 = row.get("chr1", "")
                x1   = row.get("x1", "")
                x2   = row.get("x2", "")
                chr2 = row.get("chr2", "")
                y1   = row.get("y1", "")
                y2   = row.get("y2", "")

                self.loops.append({
                    "chr1": chr1, "start1": x1, "end1": x2,
                    "chr2": chr2, "start2": y1, "end2": y2,
                    "observed": row.get("observed", ""),
                    "fdr_bl": row.get("fdrBL", ""),
                    "fdr_donut": row.get("fdrDonut", ""),
                    "fdr_h": row.get("fdrH", ""),
                    "fdr_v": row.get("fdrV", ""),
                    "centroid1": row.get("centroid1", ""),
                    "centroid2": row.get("centroid2", ""),
                    "accession": accession,
                })

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        """No dedicated nodes -- loci are referenced by edges."""
        logger.info("HiChIPdb: No dedicated nodes")
        return
        yield

    def get_edges(self):
        """Yield ChromatinInteraction edges (locus1 -> locus2)."""
        logger.info("HiChIPdb: Generating ChromatinInteraction edges...")
        count = 0
        for loop in self.loops:
            locus1 = self._locus_id(loop["chr1"], loop["start1"], loop["end1"])
            locus2 = self._locus_id(loop["chr2"], loop["start2"], loop["end2"])
            edge_id = f"hichipdb:{loop['accession']}_{count}"

            props = {
                "chr1":      self._sanitize(loop["chr1"]),
                "start1":    int(loop["start1"]) if loop["start1"] else 0,
                "end1":      int(loop["end1"])   if loop["end1"]   else 0,
                "chr2":      self._sanitize(loop["chr2"]),
                "start2":    int(loop["start2"]) if loop["start2"] else 0,
                "end2":      int(loop["end2"])   if loop["end2"]   else 0,
                "observed":  self._sanitize(loop["observed"]),
                "fdr_bl":    self._sanitize(loop["fdr_bl"]),
                "fdr_donut": self._sanitize(loop["fdr_donut"]),
                "centroid1": self._sanitize(loop["centroid1"]),
                "centroid2": self._sanitize(loop["centroid2"]),
                "accession": self._sanitize(loop["accession"]),
                "source":    "HiChIPdb",
            }
            yield (edge_id, locus1, locus2, "chromatin interaction", props)
            count += 1
        logger.info(f"HiChIPdb: Generated {count} ChromatinInteraction edges")
