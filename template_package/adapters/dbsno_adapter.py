"""
dbSNO (S-Nitrosylation Database) Adapter for BioCypher.

Loads protein S-nitrosylation site data from multiple sources and generates
ProteinHasPTM edges (ptm_type = "s-nitrosylation").

Data files parsed (in priority order):
  1. sno_sites_all_species.csv    -- curated dbSNO sites (CSV, 556 records)
  2. sno_sites_from_uniprot.csv   -- UniProt-extracted sites (CSV, 78 records)
  3. uniprot_snitrosocysteine.tsv -- UniProt Modified residue (TSV, 502 proteins)
  4. uniprot_nitrosylation_ptm.tsv -- UniProt PTM + Modified residue (TSV, 221)
  5. uniprot_sno_all_species.tsv  -- UniProt all species (TSV, 500 proteins)
  6. uniprot_snitrosocysteine_human.tsv -- human subset (TSV, 57)
  7. uniprot_sno_human.tsv        -- human subset (TSV, 58)
  8. uniprot_sno_human2.tsv       -- human subset (TSV, 149)
  9. uniprot_sno_ptm_comment.tsv  -- human PTM comment (TSV, 41)

Deduplication key: (uniprot_acc, position).
"""

import re
from pathlib import Path
from biocypher._logger import logger


class DbSNOAdapter:
    def __init__(self, data_dir="template_package/data/dbsno"):
        self.data_dir = Path(data_dir)
        self.sites = []          # list of dicts, one per individual SNO site
        self._seen_keys = set()  # dedup key = (acc, position)
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

    _MODRES_SNO_RE = re.compile(
        r'MOD_RES\s+(\d+)\s*;\s*/note="S-nitrosocysteine[^"]*"'
        r'(?:\s*;\s*/evidence="([^"]*)")?'
    )

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("dbSNO: data directory not found")
            return

        # 1. Curated CSV files (dbSNO-specific format)
        for csv_name in ["sno_sites_all_species.csv", "sno_sites_from_uniprot.csv"]:
            fpath = self.data_dir / csv_name
            if fpath.exists():
                try:
                    self._parse_dbsno_csv(fpath)
                except Exception as e:
                    logger.warning(f"dbSNO: Error reading {csv_name}: {e}")

        # 2. UniProt TSV files containing Modified residue / Lipidation columns
        tsv_files = [
            "uniprot_snitrosocysteine.tsv",
            "uniprot_nitrosylation_ptm.tsv",
            "uniprot_sno_all_species.tsv",
            "uniprot_snitrosocysteine_human.tsv",
            "uniprot_sno_human.tsv",
            "uniprot_sno_human2.tsv",
            "uniprot_sno_ptm_comment.tsv",
        ]
        for tsv_name in tsv_files:
            fpath = self.data_dir / tsv_name
            if fpath.exists():
                try:
                    self._parse_uniprot_tsv(fpath)
                except Exception as e:
                    logger.warning(f"dbSNO: Error reading {tsv_name}: {e}")

        logger.info(
            f"dbSNO: Loaded {len(self.sites)} S-nitrosylation sites "
            f"from {len({s['acc'] for s in self.sites})} proteins"
        )

    def _parse_dbsno_csv(self, fpath):
        """Parse a dbSNO curated CSV (comma-separated, with header).

        Expected columns: uniprot_id, gene_name, [organism,] protein_name,
        position, modification, evidence
        """
        with open(fpath, "r", errors="replace") as fh:
            header_line = fh.readline()
            if header_line.startswith("<"):
                return
            headers = [h.strip() for h in header_line.rstrip("\n").split(",")]
            col = {h: i for i, h in enumerate(headers)}

            acc_idx = col.get("uniprot_id")
            gene_idx = col.get("gene_name")
            pos_idx = col.get("position")
            mod_idx = col.get("modification")
            ev_idx = col.get("evidence")
            org_idx = col.get("organism")
            prot_idx = col.get("protein_name")

            if acc_idx is None or pos_idx is None:
                return

            for line in fh:
                # Handle quoted fields (evidence may contain commas)
                parts = self._split_csv_line(line.rstrip("\n"))
                if len(parts) <= max(acc_idx, pos_idx):
                    continue

                acc = parts[acc_idx].strip()
                position = parts[pos_idx].strip()
                if not acc or not position:
                    continue

                key = (acc, position)
                if key in self._seen_keys:
                    continue
                self._seen_keys.add(key)

                gene = parts[gene_idx].strip() if gene_idx is not None and gene_idx < len(parts) else ""
                organism = parts[org_idx].strip() if org_idx is not None and org_idx < len(parts) else ""
                evidence = parts[ev_idx].strip() if ev_idx is not None and ev_idx < len(parts) else ""

                # Extract PMIDs from evidence strings like "ECO:...|PubMed:12345"
                pmids = self._extract_pmids(evidence)

                self.sites.append({
                    "acc": acc,
                    "gene": gene,
                    "position": position,
                    "organism": organism,
                    "pmids": pmids,
                    "evidence": evidence,
                    "source_file": fpath.name,
                })

    @staticmethod
    def _split_csv_line(line):
        """Simple CSV split that respects double-quoted fields."""
        result = []
        current = []
        in_quotes = False
        for ch in line:
            if ch == '"':
                in_quotes = not in_quotes
            elif ch == ',' and not in_quotes:
                result.append(''.join(current))
                current = []
            else:
                current.append(ch)
        result.append(''.join(current))
        return result

    def _parse_uniprot_tsv(self, fpath):
        """Parse a UniProt-style TSV that has Entry/Gene Names columns and a
        Modified residue column containing MOD_RES annotations.  Only
        S-nitrosocysteine sites are extracted."""
        with open(fpath, "r", errors="replace") as fh:
            header_line = fh.readline()
            if header_line.startswith("<"):
                return
            headers = header_line.rstrip("\n").split("\t")
            col = {h.strip(): i for i, h in enumerate(headers)}

            acc_idx = col.get("Entry")
            gene_idx = col.get("Gene Names")
            modres_idx = col.get("Modified residue")
            organism_idx = col.get("Organism")

            if acc_idx is None or modres_idx is None:
                return

            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(acc_idx, modres_idx):
                    continue

                acc = parts[acc_idx].strip()
                if not acc:
                    continue

                modres_field = parts[modres_idx]
                gene_raw = parts[gene_idx].strip() if gene_idx is not None and gene_idx < len(parts) else ""
                gene = gene_raw.split()[0] if gene_raw else ""
                organism = parts[organism_idx].strip() if organism_idx is not None and organism_idx < len(parts) else ""

                for m in self._MODRES_SNO_RE.finditer(modres_field):
                    position = m.group(1)
                    evidence = m.group(2) or ""
                    key = (acc, position)
                    if key in self._seen_keys:
                        continue
                    self._seen_keys.add(key)

                    pmids = self._extract_pmids(evidence)

                    self.sites.append({
                        "acc": acc,
                        "gene": gene,
                        "position": position,
                        "organism": organism,
                        "pmids": pmids,
                        "evidence": evidence,
                        "source_file": fpath.name,
                    })

    @staticmethod
    def _extract_pmids(evidence_str):
        """Pull PubMed IDs from evidence strings."""
        if not evidence_str:
            return ""
        pmids = re.findall(r'PubMed:(\d+)', evidence_str)
        return ";".join(sorted(set(pmids)))

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        logger.info("dbSNO: No dedicated nodes (uses existing protein nodes)")
        return
        yield

    def get_edges(self):
        logger.info("dbSNO: Generating ProteinHasPTM edges (S-nitrosylation)...")
        count = 0
        for site in self.sites:
            edge_id = f"dbsno:{site['acc']}_C{site['position']}"
            target = site["gene"] if site["gene"] else site["acc"]
            props = {
                "site": self._sanitize(site["position"]),
                "ptm_type": "s-nitrosylation",
                "residue": "Cys",
                "score": 0,
                "enzymes": "",
                "pmids": self._sanitize(site["pmids"]),
                "source": "dbSNO",
                "organism": self._sanitize(site["organism"]),
                "evidence": self._sanitize(site["evidence"]),
            }
            yield (None, site["acc"], target, "ProteinHasPTM", props)
            count += 1
        logger.info(f"dbSNO: Generated {count} ProteinHasPTM edges")
