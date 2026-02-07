"""
SEPDB (Super-Enhancer Database / dbSUPER) Adapter for BioCypher.

Parses dbSUPER_SuperEnhancers_hg19.tsv to yield:
  - SuperEnhancer nodes  (genomic regions with tissue annotations)
  - Regulatory edges     (super enhancer -> gene)
"""

from pathlib import Path
from biocypher._logger import logger


class SEPDBAdapter:
    def __init__(self, data_dir="template_package/data/sepdb"):
        self.data_dir = Path(data_dir)
        self.enhancers = []   # list of dicts
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

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("SEPDB: data directory not found")
            return
        self._parse_tsv()
        logger.info(f"SEPDB: Loaded {len(self.enhancers)} super-enhancers")

    def _parse_tsv(self):
        """
        Parse dbSUPER_SuperEnhancers_hg19.tsv.
        Columns (note leading spaces): chrom, start, stop, se_id, gene_symbol, cell_name, rank
        """
        fpath = self.data_dir / "dbSUPER_SuperEnhancers_hg19.tsv"
        if not fpath.exists():
            logger.warning("SEPDB: dbSUPER TSV not found")
            return
        with open(fpath, "r", errors="replace") as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if header is None:
                    # strip whitespace from column names
                    header = [c.strip() for c in parts]
                    continue
                if len(parts) < 5:
                    continue
                row = dict(zip(header, [p.strip() for p in parts]))

                self.enhancers.append({
                    "chrom":       row.get("chrom", ""),
                    "start":       row.get("start", ""),
                    "stop":        row.get("stop", ""),
                    "se_id":       row.get("se_id", ""),
                    "gene_symbol": row.get("gene_symbol", ""),
                    "cell_name":   row.get("cell_name", ""),
                    "rank":        row.get("rank", ""),
                })

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        """Yield SuperEnhancer nodes."""
        logger.info("SEPDB: Generating SuperEnhancer nodes...")
        count = 0
        seen = set()
        for se in self.enhancers:
            se_id = se["se_id"]
            if not se_id or se_id in seen:
                continue
            seen.add(se_id)
            props = {
                "chrom":     self._sanitize(se["chrom"]),
                "start":     int(se["start"]) if se["start"] else 0,
                "stop":      int(se["stop"])  if se["stop"]  else 0,
                "cell_name": self._sanitize(se["cell_name"]),
                "rank":      self._sanitize(se["rank"]),
                "source":    "dbSUPER",
            }
            yield (se_id, "super enhancer", props)
            count += 1
        logger.info(f"SEPDB: Generated {count} SuperEnhancer nodes")

    def get_edges(self):
        """Yield 'super enhancer regulates gene' edges (SE -> gene symbol)."""
        logger.info("SEPDB: Generating regulatory edges...")
        count = 0
        for se in self.enhancers:
            se_id = se["se_id"]
            gene  = se["gene_symbol"]
            if not se_id or not gene:
                continue

            edge_id = f"sepdb:{se_id}_{self._sanitize(gene)}_{self._sanitize(se['cell_name'])}"
            props = {
                "cell_name": self._sanitize(se["cell_name"]),
                "rank":      self._sanitize(se["rank"]),
                "chrom":     self._sanitize(se["chrom"]),
                "start":     int(se["start"]) if se["start"] else 0,
                "stop":      int(se["stop"])  if se["stop"]  else 0,
                "source":    "dbSUPER",
            }
            yield (edge_id, se_id, gene, "super enhancer regulates gene", props)
            count += 1
        logger.info(f"SEPDB: Generated {count} regulatory edges")
