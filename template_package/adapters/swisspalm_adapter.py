"""
SwissPalm (Protein S-Palmitoylation) Adapter for BioCypher.

Loads protein palmitoylation/lipidation data from UniProt-sourced TSV files
and generates ProteinHasPTM edges with per-site granularity.

Data files parsed:
  - uniprot_palmitoylation_reviewed.tsv   (4,844 proteins, all species, reviewed)
  - uniprot_all_lipidation_human.tsv      (801 human lipidation entries)
  - uniprot_palmitoylation_human_detailed.tsv (314 human, detailed)
  - uniprot_palmitoylation.tsv            (321 records)

Each protein may carry multiple LIPID annotations; every individual site
is emitted as a separate ProteinHasPTM edge.
"""

import re
from pathlib import Path
from biocypher._logger import logger


class SwissPalmAdapter:
    def __init__(self, data_dir="template_package/data/swisspalm"):
        self.data_dir = Path(data_dir)
        self.sites = []          # list of dicts, one per individual lipid site
        self._seen_keys = set()  # dedup key = (acc, position, note)
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
        return text.strip()

    _LIPID_RE = re.compile(
        r'LIPID\s+(\d+)\s*;\s*/note="([^"]+)"'
        r'(?:\s*;\s*/evidence="([^"]*)")?'
    )

    def _parse_lipid_field(self, lipid_str):
        """Yield (position, note, evidence) tuples from a Lipidation column."""
        for m in self._LIPID_RE.finditer(lipid_str):
            yield m.group(1), m.group(2), m.group(3) or ""

    @staticmethod
    def _classify_ptm(note):
        """Return a normalised ptm_type string from the LIPID /note value."""
        nl = note.lower()
        if "palmitoyl" in nl:
            return "palmitoylation"
        if "myristoyl" in nl:
            return "myristoylation"
        if "farnesyl" in nl:
            return "prenylation"
        if "geranylgeranyl" in nl:
            return "prenylation"
        if "gpi-anchor" in nl or "gpi" in nl:
            return "GPI-anchor"
        if "diacylglycerol" in nl:
            return "lipidation"
        return "lipidation"

    @staticmethod
    def _first_gene(gene_field):
        """Return the first token of a potentially multi-gene field."""
        if not gene_field:
            return ""
        return gene_field.split()[0]

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("SwissPalm: data directory not found")
            return

        # Process files in a defined priority order so that the most
        # informative file (reviewed, all species) is loaded first.
        file_order = [
            "uniprot_palmitoylation_reviewed.tsv",
            "uniprot_all_lipidation_human.tsv",
            "uniprot_palmitoylation_human_detailed.tsv",
            "uniprot_palmitoylation.tsv",
        ]
        for fname in file_order:
            fpath = self.data_dir / fname
            if not fpath.exists():
                continue
            try:
                self._parse_uniprot_tsv(fpath)
            except Exception as e:
                logger.warning(f"SwissPalm: Error reading {fpath.name}: {e}")

        logger.info(
            f"SwissPalm: Loaded {len(self.sites)} individual lipidation "
            f"sites from {len({s['acc'] for s in self.sites})} proteins"
        )

    def _parse_uniprot_tsv(self, fpath):
        """Parse a UniProt-style TSV that has an Entry column and a
        Lipidation column containing one or more ``LIPID ...`` annotations."""
        with open(fpath, "r", errors="replace") as fh:
            header_line = fh.readline()
            if header_line.startswith("<"):
                return  # HTML error page
            headers = header_line.rstrip("\n").split("\t")

            # Identify relevant column indices
            col = {h: i for i, h in enumerate(headers)}
            acc_idx = col.get("Entry")
            gene_idx = col.get("Gene Names", col.get("Gene_name"))
            lipid_idx = col.get("Lipidation")
            organism_idx = col.get("Organism")

            if acc_idx is None or lipid_idx is None:
                return  # cannot process without these columns

            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(acc_idx, lipid_idx):
                    continue

                acc = parts[acc_idx].strip()
                if not acc:
                    continue

                gene_raw = parts[gene_idx].strip() if gene_idx is not None and gene_idx < len(parts) else ""
                gene = self._first_gene(gene_raw)
                organism = parts[organism_idx].strip() if organism_idx is not None and organism_idx < len(parts) else ""
                lipid_field = parts[lipid_idx]

                for position, note, evidence in self._parse_lipid_field(lipid_field):
                    key = (acc, position, note)
                    if key in self._seen_keys:
                        continue
                    self._seen_keys.add(key)

                    self.sites.append({
                        "acc": acc,
                        "gene": gene,
                        "position": position,
                        "note": note,
                        "evidence": evidence,
                        "organism": organism,
                        "ptm_type": self._classify_ptm(note),
                        "source_file": fpath.name,
                    })

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        logger.info("SwissPalm: No dedicated nodes (uses existing protein nodes)")
        return
        yield

    def get_edges(self):
        logger.info("SwissPalm: Generating ProteinHasPTM edges...")
        count = 0
        for site in self.sites:
            edge_id = f"swisspalm:{site['acc']}_LIPID{site['position']}"
            target = site["gene"] if site["gene"] else site["acc"]
            props = {
                "site": self._sanitize(site["position"]),
                "ptm_type": site["ptm_type"],
                "residue": self._sanitize(site["note"]),
                "score": 0,
                "enzymes": "",
                "pmids": "",
                "source": "SwissPalm",
                "organism": self._sanitize(site["organism"]),
                "evidence": self._sanitize(site["evidence"]),
            }
            yield (None, site["acc"], target, "ProteinHasPTM", props)
            count += 1
        logger.info(f"SwissPalm: Generated {count} ProteinHasPTM edges")
