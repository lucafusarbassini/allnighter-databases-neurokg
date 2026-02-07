"""
UbiNet (E3/DUB-Substrate Interaction) Adapter for BioCypher.

Loads literature-curated E3 ubiquitin-ligase and deubiquitinase (DUB)
substrate interaction data and generates ProteinHasPTM edges
(ubiquitination / deubiquitination).

Data files parsed:
  - literature.E3.txt   (4,068 E3-substrate interactions)
  - literature.DUB.txt  (967 DUB-substrate interactions)

Columns (tab-separated, identical layout for E3 and DUB files):
  NUMBER | SwissProt ID (E3/DUB) | SwissProt ID (Substrate) |
  SwissProt AC (E3/DUB) | SwissProt AC (Substrate) |
  Gene Symbol (E3/DUB) | Gene Symbol (Substrate) |
  SOURCE | SOURCEID | SENTENCE | E3TYPE/DUBTYPE | COUNT | type | species
"""

from pathlib import Path
from biocypher._logger import logger


class UbiNetAdapter:
    def __init__(self, data_dir="template_package/data/ubinet"):
        self.data_dir = Path(data_dir)
        self.interactions = []   # list of dicts
        self._seen_keys = set()  # dedup key = (enzyme_ac, substrate_ac, interaction_type)
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
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        # Strip HTML tags that appear in some SENTENCE fields
        import re
        text = re.sub(r'<[^>]+>', '', text)
        return text.strip()

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("UbiNet: data directory not found")
            return

        for fname, interaction_type in [
            ("literature.E3.txt", "E3-substrate"),
            ("literature.DUB.txt", "DUB-substrate"),
        ]:
            fpath = self.data_dir / fname
            if not fpath.exists():
                logger.warning(f"UbiNet: {fname} not found")
                continue
            try:
                self._parse_interaction_file(fpath, interaction_type)
            except Exception as e:
                logger.warning(f"UbiNet: Error reading {fname}: {e}")

        logger.info(f"UbiNet: Loaded {len(self.interactions)} interactions")

    def _parse_interaction_file(self, fpath, interaction_type):
        """Parse a UbiNet literature interaction file (tab-separated)."""
        with open(fpath, "r", errors="replace") as fh:
            header_line = fh.readline()
            if header_line.startswith("<"):
                return
            headers = header_line.rstrip("\n").split("\t")
            col = {h.strip(): i for i, h in enumerate(headers)}

            # Identify column indices -- use the actual header names
            # E3 file: "SwissProt AC (E3)", DUB file: "SwissProt AC (DUB)"
            enzyme_ac_idx = None
            enzyme_gene_idx = None
            enzyme_type_idx = None
            for h, i in col.items():
                if h.startswith("SwissProt AC (E3") or h.startswith("SwissProt AC (DUB"):
                    enzyme_ac_idx = i
                if h.startswith("Gene Symbol (E3") or h.startswith("Gene Symbol (DUB"):
                    enzyme_gene_idx = i
                if h == "E3TYPE" or h == "DUBTYPE":
                    enzyme_type_idx = i

            substrate_ac_idx = col.get("SwissProt AC (Substrate)")
            substrate_gene_idx = col.get("Gene Symbol (Substrate)")
            source_idx = col.get("SOURCE")
            sourceid_idx = col.get("SOURCEID")
            species_idx = col.get("species")
            type_idx = col.get("type")

            if enzyme_ac_idx is None or substrate_ac_idx is None:
                logger.warning(f"UbiNet: Could not find required columns in {fpath.name}")
                return

            for line in fh:
                parts = line.rstrip("\n").split("\t")
                max_needed = max(enzyme_ac_idx, substrate_ac_idx)
                if len(parts) <= max_needed:
                    continue

                enzyme_ac = parts[enzyme_ac_idx].strip()
                substrate_ac = parts[substrate_ac_idx].strip()
                if not enzyme_ac or not substrate_ac:
                    continue

                key = (enzyme_ac, substrate_ac, interaction_type)
                if key in self._seen_keys:
                    continue
                self._seen_keys.add(key)

                enzyme_gene = parts[enzyme_gene_idx].strip() if enzyme_gene_idx is not None and enzyme_gene_idx < len(parts) else ""
                substrate_gene = parts[substrate_gene_idx].strip() if substrate_gene_idx is not None and substrate_gene_idx < len(parts) else ""
                etype = parts[enzyme_type_idx].strip() if enzyme_type_idx is not None and enzyme_type_idx < len(parts) else ""
                source = parts[source_idx].strip() if source_idx is not None and source_idx < len(parts) else ""
                sourceid = parts[sourceid_idx].strip() if sourceid_idx is not None and sourceid_idx < len(parts) else ""
                species = parts[species_idx].strip() if species_idx is not None and species_idx < len(parts) else ""
                int_subtype = parts[type_idx].strip() if type_idx is not None and type_idx < len(parts) else ""

                self.interactions.append({
                    "enzyme_ac": enzyme_ac,
                    "enzyme_gene": enzyme_gene,
                    "substrate_ac": substrate_ac,
                    "substrate_gene": substrate_gene,
                    "interaction_type": interaction_type,
                    "enzyme_class": etype,
                    "source": source,
                    "pmid": sourceid,
                    "species": species,
                    "subtype": int_subtype,
                })

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        logger.info("UbiNet: No dedicated nodes (uses existing protein nodes)")
        return
        yield

    def get_edges(self):
        logger.info("UbiNet: Generating ProteinHasPTM edges (ubiquitination)...")
        count = 0
        for ix in self.interactions:
            edge_id = f"ubinet:{ix['enzyme_ac']}_{ix['substrate_ac']}_{ix['interaction_type']}"
            ptm_type = "ubiquitination" if ix["interaction_type"] == "E3-substrate" else "deubiquitination"
            props = {
                "site": "",
                "ptm_type": ptm_type,
                "residue": "",
                "score": 0,
                "enzymes": self._sanitize(ix["enzyme_gene"]),
                "enzyme_uniprot": ix["enzyme_ac"],
                "enzyme_class": self._sanitize(ix["enzyme_class"]),
                "interaction_type": ix["interaction_type"],
                "pmids": self._sanitize(ix["pmid"]),
                "source": "UbiNet",
                "species": self._sanitize(ix["species"]),
            }
            # source=enzyme protein, target=substrate protein
            yield (
                None,
                ix["substrate_ac"],
                ix["substrate_gene"] if ix["substrate_gene"] else ix["substrate_ac"],
                "ProteinHasPTM",
                props,
            )
            count += 1
        logger.info(f"UbiNet: Generated {count} ProteinHasPTM edges")
