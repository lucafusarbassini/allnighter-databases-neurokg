"""
Proteoform Atlas Adapter for BioCypher.
Parses UniProt human proteoform data from three TSV files:
  - uniprot_human_variants.tsv        (splice variants / VAR_SEQ)
  - uniprot_human_ptm_processing.tsv  (chain, signal peptide, transit peptide, etc.)
  - uniprot_human_isoforms_full.tsv   (alternative products / isoform descriptions)

Yields:
  - Proteoform variant NODES  (splice variants, chains, signal peptides, etc.)
  - Protein-to-proteoform EDGES linking UniProt accession to each variant/feature
"""

import csv
import hashlib
import re
from pathlib import Path
from biocypher._logger import logger


class ProteoformAtlasAdapter:
    def __init__(self, data_dir="template_package/data/proteoform_atlas"):
        self.data_dir = Path(data_dir)
        self.variants = []        # parsed from variants TSV
        self.features = []        # parsed from ptm_processing TSV
        self.isoforms = []        # parsed from isoforms TSV
        self._load_data()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    @staticmethod
    def _make_id(*parts):
        raw = "|".join(str(p) for p in parts)
        return hashlib.md5(raw.encode()).hexdigest()[:12]

    # ------------------------------------------------------------------
    # Feature-string parser (shared between variants and processing)
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_position(pos_str):
        """
        Parse a UniProt position string like '1..515', '40', '1..?', '?..100'.
        Returns dict with position, position_start, position_end.
        """
        result = {"position": pos_str}
        if ".." in pos_str:
            start_s, end_s = pos_str.split("..", 1)
            start_s = start_s.strip()
            end_s = end_s.strip()
            if start_s.isdigit():
                result["position_start"] = int(start_s)
            if end_s.isdigit():
                result["position_end"] = int(end_s)
        else:
            pos_s = pos_str.strip()
            if pos_s.isdigit():
                result["position_start"] = int(pos_s)
                result["position_end"] = int(pos_s)
        return result

    @staticmethod
    def _parse_features(raw_text, feature_prefix):
        """
        Parse UniProt feature strings like:
          'VAR_SEQ 40; /note="Missing (in isoform 2)"; /id="VSP_062288"'
        into a list of dicts.  Multiple features in one cell are separated by
        '; ' followed by an uppercase keyword.
        """
        if not raw_text or not raw_text.strip():
            return []
        results = []
        # Split on boundaries between features
        parts = re.split(r';\s+(?=[A-Z_]+ )', raw_text)
        for part in parts:
            part = part.strip().rstrip(';').strip()
            if not part:
                continue
            feat = {"raw": part}
            # Extract feature type and position (supports '?', e.g. '1..?')
            m = re.match(r"([A-Z_]+)\s+([\d?]+\.\.[\d?]+|[\d]+)", part)
            if m:
                feat["type"] = m.group(1)
                pos_info = ProteoformAtlasAdapter._parse_position(m.group(2))
                feat.update(pos_info)
            else:
                feat["type"] = feature_prefix
            # Extract /note
            note_m = re.search(r'/note="([^"]*)"', part)
            if note_m:
                feat["description"] = note_m.group(1)
            # Extract /id
            id_m = re.search(r'/id="([^"]*)"', part)
            if id_m:
                feat["feature_id"] = id_m.group(1)
            # Extract /evidence
            ev_m = re.search(r'/evidence="([^"]*)"', part)
            if ev_m:
                feat["evidence"] = ev_m.group(1)

            results.append(feat)
        return results

    # ------------------------------------------------------------------
    # data loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("ProteoformAtlasAdapter: data directory not found")
            return

        self._load_variants()
        self._load_ptm_processing()
        self._load_isoforms()

        total = len(self.variants) + len(self.features) + len(self.isoforms)
        logger.info(
            f"ProteoformAtlasAdapter: Loaded {len(self.variants)} splice variants, "
            f"{len(self.features)} processing features, "
            f"{len(self.isoforms)} isoform records ({total} total)"
        )

    def _load_tsv(self, filename):
        """Generic TSV loader; returns list of dicts."""
        path = self.data_dir / filename
        if not path.exists():
            logger.warning(f"ProteoformAtlasAdapter: {filename} not found")
            return []
        rows = []
        with open(path, "r", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(row)
        return rows

    def _load_variants(self):
        """Parse splice variant records from the variants TSV."""
        rows = self._load_tsv("uniprot_human_variants.tsv")
        for row in rows:
            entry = row.get("Entry", "").strip()
            gene = row.get("Gene Names", "").strip()
            protein = row.get("Protein names", "").strip()
            alt_seq = row.get("Alternative sequence", "").strip()
            if not entry or not alt_seq:
                continue
            features = self._parse_features(alt_seq, "VAR_SEQ")
            for feat in features:
                feat["entry"] = entry
                feat["gene"] = gene
                feat["protein_name"] = protein
                feat["variant_type"] = "splice_variant"
                self.variants.append(feat)

    def _load_ptm_processing(self):
        """Parse protein processing features (chain, signal, transit, etc.)."""
        rows = self._load_tsv("uniprot_human_ptm_processing.tsv")
        feature_cols = [
            "Chain", "Initiator methionine", "Signal peptide",
            "Transit peptide", "Propeptide", "Peptide",
        ]
        for row in rows:
            entry = row.get("Entry", "").strip()
            gene = row.get("Gene Names", "").strip()
            protein = row.get("Protein names", "").strip()
            if not entry:
                continue
            for col in feature_cols:
                cell = row.get(col, "").strip()
                if not cell:
                    continue
                features = self._parse_features(cell, col.upper().replace(" ", "_"))
                for feat in features:
                    feat["entry"] = entry
                    feat["gene"] = gene
                    feat["protein_name"] = protein
                    feat["variant_type"] = col.lower().replace(" ", "_")
                    self.features.append(feat)

    def _load_isoforms(self):
        """Parse isoform records."""
        rows = self._load_tsv("uniprot_human_isoforms_full.tsv")
        for row in rows:
            entry = row.get("Entry", "").strip()
            gene = row.get("Gene Names", "").strip()
            protein = row.get("Protein names", "").strip()
            alt_prod = row.get("Alternative products (isoforms)", "").strip()
            if not entry or not alt_prod:
                continue
            # Extract named isoforms
            isoform_names = re.findall(r"Name=([^;]+)", alt_prod)
            isoform_ids = re.findall(r"IsoId=([^;]+)", alt_prod)
            n_isoforms_m = re.search(r"Named isoforms=(\d+)", alt_prod)
            n_isoforms = int(n_isoforms_m.group(1)) if n_isoforms_m else len(isoform_names)
            events = re.findall(r"Event=([^;]+)", alt_prod)

            rec = {
                "entry": entry,
                "gene": gene,
                "protein_name": protein,
                "variant_type": "isoform",
                "n_isoforms": n_isoforms,
                "isoform_names": "; ".join(isoform_names),
                "isoform_ids": "; ".join(isoform_ids),
                "events": "; ".join(events),
                "description": alt_prod[:500],
            }
            self.isoforms.append(rec)

    # ------------------------------------------------------------------
    # BioCypher interface: NODES
    # ------------------------------------------------------------------
    def get_nodes(self):
        """
        Yield proteoform variant/feature nodes.
        Each node represents a distinct protein variant or processing feature.
        """
        count = 0

        # 1. Splice variant nodes
        for feat in self.variants:
            fid = feat.get("feature_id", "")
            node_id = f"splice_variant:{fid}" if fid else \
                f"splice_variant:{self._make_id(feat['entry'], feat.get('position',''), feat.get('description',''))}"
            props = {
                "uniprot_entry": feat["entry"],
                "gene": self._sanitize(feat.get("gene", "")),
                "protein_name": self._sanitize(feat.get("protein_name", ""))[:200],
                "variant_type": "splice_variant",
                "position": feat.get("position", ""),
                "description": self._sanitize(feat.get("description", "")),
                "feature_id": fid,
                "evidence": feat.get("evidence", ""),
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield node_id, "proteoform_variant", props

        # 2. Processing feature nodes (chain, signal peptide, etc.)
        for feat in self.features:
            fid = feat.get("feature_id", "")
            node_id = f"processing:{fid}" if fid else \
                f"processing:{self._make_id(feat['entry'], feat.get('type',''), feat.get('position',''))}"
            props = {
                "uniprot_entry": feat["entry"],
                "gene": self._sanitize(feat.get("gene", "")),
                "protein_name": self._sanitize(feat.get("protein_name", ""))[:200],
                "variant_type": feat.get("variant_type", ""),
                "feature_type": feat.get("type", ""),
                "position": feat.get("position", ""),
                "description": self._sanitize(feat.get("description", "")),
                "feature_id": fid,
                "evidence": feat.get("evidence", ""),
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield node_id, "proteoform_feature", props

        # 3. Isoform nodes
        for rec in self.isoforms:
            node_id = f"isoform_set:{rec['entry']}"
            props = {
                "uniprot_entry": rec["entry"],
                "gene": self._sanitize(rec.get("gene", "")),
                "protein_name": self._sanitize(rec.get("protein_name", ""))[:200],
                "variant_type": "isoform",
                "n_isoforms": rec.get("n_isoforms", 0),
                "isoform_names": self._sanitize(rec.get("isoform_names", "")),
                "isoform_ids": self._sanitize(rec.get("isoform_ids", "")),
                "events": self._sanitize(rec.get("events", "")),
                "description": self._sanitize(rec.get("description", ""))[:300],
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield node_id, "proteoform_isoform", props

        logger.info(f"ProteoformAtlasAdapter: Yielded {count} proteoform nodes")

    # ------------------------------------------------------------------
    # BioCypher interface: EDGES
    # ------------------------------------------------------------------
    def get_edges(self):
        """
        Yield edges linking UniProt protein entries to their proteoform variants.
        protein --[protein_has_variant]--> proteoform node
        """
        count = 0

        # 1. Splice variant edges
        for feat in self.variants:
            fid = feat.get("feature_id", "")
            target_id = f"splice_variant:{fid}" if fid else \
                f"splice_variant:{self._make_id(feat['entry'], feat.get('position',''), feat.get('description',''))}"
            props = {
                "variant_type": "splice_variant",
                "position": feat.get("position", ""),
                "description": self._sanitize(feat.get("description", "")),
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield (
                None,
                feat["entry"],
                target_id,
                "protein_has_variant",
                props,
            )

        # 2. Processing feature edges
        for feat in self.features:
            fid = feat.get("feature_id", "")
            target_id = f"processing:{fid}" if fid else \
                f"processing:{self._make_id(feat['entry'], feat.get('type',''), feat.get('position',''))}"
            props = {
                "variant_type": feat.get("variant_type", ""),
                "position": feat.get("position", ""),
                "description": self._sanitize(feat.get("description", "")),
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield (
                None,
                feat["entry"],
                target_id,
                "protein_has_feature",
                props,
            )

        # 3. Isoform edges
        for rec in self.isoforms:
            target_id = f"isoform_set:{rec['entry']}"
            props = {
                "variant_type": "isoform",
                "n_isoforms": rec.get("n_isoforms", 0),
                "events": self._sanitize(rec.get("events", "")),
                "source": "UniProt_proteoform_atlas",
            }
            count += 1
            yield (
                None,
                rec["entry"],
                target_id,
                "protein_has_isoforms",
                props,
            )

        logger.info(f"ProteoformAtlasAdapter: Yielded {count} proteoform edges")
