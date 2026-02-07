"""
Cistrome Data Browser (ChIP-Seq) Adapter for BioCypher.

Parses cistrome_dc2_samples.tsv, cistrome_dc2_samples.json, and
cistrome_dc2_extended_samples.json to yield ChIPSeqSample nodes
with transcription factor, cell type, tissue, and disease metadata.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class CistromeAdapter:
    def __init__(self, data_dir="template_package/data/cistrome"):
        self.data_dir = Path(data_dir)
        self.samples = {}   # id -> sample dict (deduplicated)
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
            logger.warning("Cistrome: data directory not found")
            return
        self._load_tsv()
        self._load_json(self.data_dir / "cistrome_dc2_samples.json")
        self._load_json(self.data_dir / "cistrome_dc2_extended_samples.json")
        logger.info(f"Cistrome: Loaded {len(self.samples)} unique samples")

    def _load_tsv(self):
        """Parse cistrome_dc2_samples.tsv."""
        fpath = self.data_dir / "cistrome_dc2_samples.tsv"
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
                if len(parts) < 2:
                    continue
                row = dict(zip(header, parts))
                sid = row.get("id", "")
                if not sid:
                    continue
                self.samples[str(sid)] = {
                    "id":        sid,
                    "factor":    row.get("factor", ""),
                    "cell_line": row.get("cell_line", ""),
                    "cell_type": row.get("cell_type", ""),
                    "tissue":    row.get("tissue", ""),
                    "species":   row.get("species", ""),
                    "disease":   row.get("disease", ""),
                    "geo_id":    row.get("geo_id", ""),
                    "pmid":      row.get("pmid", ""),
                    "journal":   row.get("journal", ""),
                    "lab":       row.get("lab", ""),
                    "peaks":     row.get("peaks", ""),
                    "frip":      row.get("frip", ""),
                }

    def _load_json(self, fpath):
        """Parse a JSON sample list -- merges with existing entries."""
        if not fpath.exists():
            return
        with open(fpath, "r", errors="replace") as fh:
            data = json.load(fh)
        if not isinstance(data, list):
            return
        for rec in data:
            sid = str(rec.get("id", ""))
            if not sid:
                continue
            if sid in self.samples:
                # merge any extra keys from JSON that TSV may lack
                for k, v in rec.items():
                    key = str(k)
                    if key not in self.samples[sid] or not self.samples[sid][key]:
                        self.samples[sid][key] = v
            else:
                self.samples[sid] = {
                    "id":        sid,
                    "factor":    rec.get("factor", ""),
                    "cell_line": rec.get("cell_line", ""),
                    "cell_type": rec.get("cell_type", ""),
                    "tissue":    rec.get("tissue", ""),
                    "species":   rec.get("species", ""),
                    "disease":   rec.get("disease", ""),
                    "geo_id":    rec.get("geo_id", ""),
                    "pmid":      rec.get("pmid", ""),
                    "journal":   rec.get("journal", ""),
                    "lab":       rec.get("lab", ""),
                    "peaks":     rec.get("peaks", ""),
                    "frip":      rec.get("frip", ""),
                }

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        """Yield ChIPSeqSample nodes."""
        logger.info("Cistrome: Generating ChIPSeqSample nodes...")
        count = 0
        for sid, s in sorted(self.samples.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 0):
            node_id = f"cistrome:{sid}"
            props = {
                "factor":    self._sanitize(s.get("factor", "")),
                "cell_line": self._sanitize(s.get("cell_line", "")),
                "cell_type": self._sanitize(s.get("cell_type", "")),
                "tissue":    self._sanitize(s.get("tissue", "")),
                "species":   self._sanitize(s.get("species", "")),
                "disease":   self._sanitize(s.get("disease", "")),
                "geo_id":    self._sanitize(s.get("geo_id", "")),
                "pmid":      self._sanitize(str(s.get("pmid", ""))),
                "journal":   self._sanitize(s.get("journal", "")),
                "lab":       self._sanitize(s.get("lab", "")),
                "peaks":     int(s["peaks"]) if str(s.get("peaks", "")).isdigit() else 0,
                "frip":      self._sanitize(str(s.get("frip", ""))),
                "source":    "Cistrome",
            }
            yield (node_id, "chip-seq sample", props)
            count += 1
        logger.info(f"Cistrome: Generated {count} ChIPSeqSample nodes")

    def get_edges(self):
        """
        Yield 'tf chip-seq binding' edges: TF factor -> sample.
        This links each transcription factor to the sample that profiled it.
        """
        logger.info("Cistrome: Generating TF binding edges...")
        count = 0
        for sid, s in sorted(self.samples.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 0):
            factor = s.get("factor", "")
            if not factor:
                continue
            sample_node = f"cistrome:{sid}"
            edge_id = f"cistrome:edge_{sid}"
            props = {
                "cell_line": self._sanitize(s.get("cell_line", "")),
                "cell_type": self._sanitize(s.get("cell_type", "")),
                "tissue":    self._sanitize(s.get("tissue", "")),
                "disease":   self._sanitize(s.get("disease", "")),
                "peaks":     int(s["peaks"]) if str(s.get("peaks", "")).isdigit() else 0,
                "source":    "Cistrome",
            }
            yield (edge_id, self._sanitize(factor), sample_node, "tf chip-seq binding", props)
            count += 1
        logger.info(f"Cistrome: Generated {count} TF binding edges")
