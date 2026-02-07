"""
iUUCD (Ubiquitin and Ubiquitin-like Conjugation Database) Adapter for BioCypher.

Parses iUUCD enzyme data and UniProt ubiquitin conjugation pathway annotations.

Yields ProteinHasPTM edges representing ubiquitination-related modifications.
"""

from pathlib import Path
from biocypher._logger import logger


class IUUCDAdapter:
    def __init__(self, data_dir="template_package/data/iuucd"):
        self.data_dir = Path(data_dir)
        self.entries = []
        self._seen_keys = set()
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _classify_enzyme(self, family_str):
        """Classify enzyme into E1/E2/E3/DUB from the iUUCD family field."""
        if not family_str:
            return "unknown"
        fl = family_str.lower()
        if 'e1' in fl or 'activating' in fl:
            return "E1"
        elif 'e2' in fl or 'conjugating' in fl:
            return "E2"
        elif 'e3' in fl or 'ligase' in fl:
            return "E3"
        elif 'dub' in fl or 'deubiquitin' in fl or 'hydrolase' in fl or 'usp' in fl:
            return "DUB"
        return "unknown"

    def _classify_enzyme_from_uniprot(self, protein_name, gene_names, go_mf):
        """Classify enzyme class from UniProt fields."""
        combined = (protein_name + " " + gene_names + " " + go_mf).lower()
        if 'ubiquitin-activating' in combined or 'uba' in gene_names.lower().split()[0:1]:
            return "E1"
        if 'ubiquitin-conjugating' in combined or 'ube2' in combined:
            return "E2"
        if 'ubiquitin-protein ligase' in combined or 'ubiquitin ligase' in combined:
            return "E3"
        if 'deubiquitin' in combined or 'ubiquitin carboxyl-terminal hydrolase' in combined:
            return "DUB"
        # Check gene name prefixes
        first_gene = gene_names.split()[0].upper() if gene_names.strip() else ""
        if first_gene.startswith('UBA') or first_gene.startswith('UBE1'):
            return "E1"
        if first_gene.startswith('UBE2') or first_gene.startswith('UBC'):
            return "E2"
        if first_gene.startswith('UBE3') or first_gene.startswith('RNF') or first_gene.startswith('TRIM'):
            return "E3"
        if first_gene.startswith('USP') or first_gene.startswith('OTUB') or first_gene.startswith('UCHL'):
            return "DUB"
        return "unknown"

    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("iUUCD: data directory not found")
            return

        # 1. Parse iUUCD native data
        iuucd_file = self.data_dir / "iuucd_human.tsv"
        if iuucd_file.exists():
            self._parse_iuucd_native(iuucd_file)
        else:
            logger.warning("iUUCD: iuucd_human.tsv not found")

        # 2. Parse UniProt ubiquitin conjugation data
        uniprot_file = self.data_dir / "uniprot_ubiquitin_conjugation.tsv"
        if uniprot_file.exists():
            self._parse_uniprot_ubiquitin(uniprot_file)
        else:
            logger.warning("iUUCD: uniprot_ubiquitin_conjugation.tsv not found")

        logger.info(f"iUUCD: Loaded {len(self.entries)} unique entries")

    def _parse_iuucd_native(self, fpath):
        """Parse iUUCD native TSV (status, iuucd_id, view_url, gene_alias_name,
        ensembl_gene_id, species, family)."""
        with open(fpath, 'r', errors='replace') as fh:
            header_line = fh.readline().strip()
            if not header_line:
                return
            headers = header_line.split('\t')

            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                row = dict(zip(headers, parts))

                iuucd_id = row.get('iuucd_id', '').strip()
                gene_aliases = row.get('gene_alias_name', '').strip()
                # Take first gene alias as primary
                gene_name = gene_aliases.split(';')[0].strip() if gene_aliases else ""
                ensembl_id = row.get('ensembl_gene_id', '').strip()
                # Strip version suffix from Ensembl ID
                if '.' in ensembl_id:
                    ensembl_id = ensembl_id.split('.')[0]
                family = row.get('family', '').strip()
                status = row.get('status', '').strip()

                if not iuucd_id or not gene_name:
                    continue

                dedup_key = ("iuucd", iuucd_id)
                if dedup_key in self._seen_keys:
                    continue
                self._seen_keys.add(dedup_key)

                enzyme_class = self._classify_enzyme(family)

                self.entries.append({
                    'id': iuucd_id,
                    'gene_name': gene_name,
                    'all_gene_aliases': gene_aliases,
                    'ensembl_id': ensembl_id,
                    'protein_name': '',
                    'enzyme_class': enzyme_class,
                    'family': family,
                    'ptm_type': 'ubiquitination',
                    'status': status,
                    'source_file': 'iuucd_human',
                })

    def _parse_uniprot_ubiquitin(self, fpath):
        """Parse UniProt ubiquitin conjugation pathway TSV."""
        with open(fpath, 'r', errors='replace') as fh:
            header_line = fh.readline().strip()
            if not header_line or header_line.startswith('<'):
                return
            headers = header_line.split('\t')

            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                row = dict(zip(headers, parts))

                uniprot_id = row.get('Entry', '').strip()
                if not uniprot_id:
                    continue

                gene_names_raw = row.get('Gene Names', '').strip()
                gene_name = gene_names_raw.split()[0] if gene_names_raw else ""
                protein_name = row.get('Protein names', '').strip()
                go_mf = row.get('Gene Ontology (molecular function)', '').strip()
                go_bp = row.get('Gene Ontology (biological process)', '').strip()
                domain = row.get('Domain [FT]', '').strip()

                dedup_key = ("uniprot", uniprot_id)
                if dedup_key in self._seen_keys:
                    continue
                self._seen_keys.add(dedup_key)

                enzyme_class = self._classify_enzyme_from_uniprot(
                    protein_name, gene_names_raw, go_mf
                )

                self.entries.append({
                    'id': uniprot_id,
                    'gene_name': gene_name,
                    'all_gene_aliases': gene_names_raw,
                    'ensembl_id': '',
                    'protein_name': protein_name,
                    'enzyme_class': enzyme_class,
                    'family': '',
                    'ptm_type': 'ubiquitin conjugation',
                    'status': '',
                    'go_molecular_function': go_mf,
                    'go_biological_process': go_bp,
                    'domain': domain,
                    'source_file': 'uniprot_ubiquitin_conjugation',
                })

    def get_nodes(self):
        logger.info("iUUCD: No dedicated nodes (proteins referenced by ID)")
        return
        yield

    def get_edges(self):
        logger.info("iUUCD: Generating ProteinHasPTM edges (ubiquitination)...")
        count = 0
        for entry in self.entries:
            edge_id = f"iuucd:{entry['id']}"
            source_id = entry['id']
            target_id = f"ptm:ubiquitination:{entry['gene_name']}"

            props = {
                'ptm_type': self._sanitize(entry['ptm_type']),
                'enzyme_class': self._sanitize(entry['enzyme_class']),
                'gene_name': self._sanitize(entry['gene_name']),
                'protein_name': self._sanitize(entry.get('protein_name', '')),
                'family': self._sanitize(entry.get('family', '')),
                'ensembl_id': self._sanitize(entry.get('ensembl_id', '')),
                'source': 'iUUCD',
                'source_file': entry.get('source_file', ''),
            }
            yield (edge_id, source_id, target_id,
                   "ProteinHasPTM", props)
            count += 1
        logger.info(f"iUUCD: Generated {count} ProteinHasPTM edges")
