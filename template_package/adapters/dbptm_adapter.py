"""
dbPTM Adapter for BioCypher.

Loads dbPTM post-translational modification data and generates:
- ProteinHasPTM edges (protein -> PTM site with type, residue, and evidence)

dbPTM is a comprehensive resource for experimentally verified and
computationally predicted post-translational modifications including
phosphorylation, ubiquitination, acetylation, methylation, SUMOylation,
glycosylation, and many other modification types.

Expected data directory: template_package/data/dbptm/
Expected files (any of):
  - PTM_of_Human.txt
  - Human_PTM.tsv or Human_PTM.tsv.gz
  - dbPTM_Human.tsv or dbPTM_Human.tsv.gz

Also supports UniProt MOD_RES flat annotation format as a fallback, where
the 'Modified residue' column contains concatenated MOD_RES annotation
strings (e.g. from a UniProt TSV export).

dbPTM typical columns (tab-separated):
  UniProtKB_AC, Position, PTM_type, Residue, Modified_residue, Source, PMIDs

Reference:
  Huang et al., dbPTM in 2019: exploring disease association and
  cross-talk of post-translational modifications.
  Nucleic Acids Res. 2019 Jan;47(D1):D298-D308.
"""

import gzip
import re
from pathlib import Path
from biocypher._logger import logger


class DbPTMAdapter:
    def __init__(self, data_dir="template_package/data/dbptm"):
        self.data_dir = Path(data_dir)
        self.ptms = []
        self._load_data()

    def _sanitize(self, text):
        """Clean text for safe use in BioCypher properties."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """
        Locate and load dbPTM data files.

        Searches for .tsv.gz, .txt.gz, .tsv, and .txt files. Each candidate
        is tested: HTML pages and corrupt gzip files are skipped. The file
        format is auto-detected (standard dbPTM columnar vs. UniProt MOD_RES
        flat annotation).
        """
        if not self.data_dir.exists():
            logger.warning("dbPTM: data directory not found")
            return

        candidates = (
            list(self.data_dir.glob("*.tsv.gz"))
            + list(self.data_dir.glob("*.txt.gz"))
            + list(self.data_dir.glob("*.tsv"))
            + list(self.data_dir.glob("*.txt"))
        )

        if not candidates:
            logger.warning("dbPTM: no data files found in data directory")
            return

        for fpath in candidates:
            try:
                if str(fpath).endswith('.gz'):
                    with gzip.open(fpath, 'rt', errors='replace') as f:
                        first = f.readline()
                        if first.lstrip().startswith('<'):
                            logger.warning(
                                f"dbPTM: {fpath.name} appears to be HTML "
                                "or corrupt, skipping"
                            )
                            continue
                        f.seek(0)
                        self._parse_auto(f, fpath.name)
                else:
                    with open(fpath, 'r', errors='replace') as f:
                        first = f.readline()
                        if first.lstrip().startswith('<'):
                            logger.warning(
                                f"dbPTM: {fpath.name} appears to be HTML, "
                                "skipping"
                            )
                            continue
                        f.seek(0)
                        self._parse_auto(f, fpath.name)
            except (gzip.BadGzipFile, OSError) as e:
                logger.warning(
                    f"dbPTM: Cannot read {fpath.name} (corrupt?): {e}"
                )
            except Exception as e:
                logger.warning(f"dbPTM: Error reading {fpath.name}: {e}")

        logger.info(f"dbPTM: Loaded {len(self.ptms)} PTM sites total")

    def _parse_auto(self, fh, filename=""):
        """
        Auto-detect the file format and dispatch to the correct parser.

        If the header contains 'Modified residue' and 'Entry' (UniProt TSV
        export format), the file is parsed as UniProt MOD_RES annotations.
        Otherwise it is parsed as standard dbPTM columnar format.
        """
        header_line = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            header_line = line
            break

        if not header_line:
            return

        header_lower = header_line.lower()

        # Detect UniProt MOD_RES format:
        # columns include 'Entry' and 'Modified residue', and data contains
        # MOD_RES annotation strings rather than simple position/type columns
        if 'modified residue' in header_lower and 'entry' in header_lower:
            self._parse_uniprot_modres(header_line, fh, filename)
        else:
            self._parse_dbptm_tabular(header_line, fh, filename)

    def _parse_dbptm_tabular(self, header_line, fh, filename=""):
        """
        Parse standard dbPTM columnar format.

        Expected columns (tab-separated):
          UniProtKB_AC, Position, PTM_type, Residue, Modified_residue,
          Source, PMIDs
        """
        header = header_line.split('\t')
        count = 0

        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('<!DOCTYPE') or line.startswith('<html'):
                logger.warning(
                    f"dbPTM: {filename} contains HTML content, stopping parse"
                )
                return

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            row = dict(zip(header, parts))

            acc = (
                row.get('UniProtKB_AC', '')
                or row.get('UniProt_AC', '')
                or row.get('uniprot_ac', '')
                or row.get('uniprot_id', '')
                or row.get('ACC', '')
                or row.get('Protein', '')
                or row.get('accession', '')
                or row.get('Entry', '')
            ).strip()

            if not acc:
                continue

            position = (
                row.get('Position', '')
                or row.get('position', '')
                or row.get('Site', '')
            ).strip()

            ptm_type = (
                row.get('PTM_type', '')
                or row.get('PTM type', '')
                or row.get('ptm_type', '')
                or row.get('Modification', '')
                or row.get('Type', '')
            ).strip()

            residue = (
                row.get('Residue', '')
                or row.get('residue', '')
                or row.get('AA', '')
            ).strip()

            pmids = (
                row.get('PMIDs', '')
                or row.get('PMID', '')
                or row.get('pmids', '')
                or row.get('PubMed', '')
                or row.get('Reference', '')
                or row.get('References', '')
            ).strip()

            ptm_type_clean = self._normalise_ptm_type(ptm_type)

            self.ptms.append({
                'acc': acc,
                'position': position,
                'ptm_type': ptm_type_clean,
                'residue': residue,
                'pmids': pmids,
            })
            count += 1

        if count > 0:
            logger.info(
                f"dbPTM: Parsed {count} records from {filename} "
                "(tabular format)"
            )

    # Regex for parsing UniProt MOD_RES annotation strings, e.g.:
    #   MOD_RES 44; /note="N6-acetyllysine"; /evidence="ECO:..."
    _MODRES_PATTERN = re.compile(
        r'MOD_RES\s+(\d+);\s*/note="([^"]*)"'
        r'(?:;\s*/evidence="([^"]*)")?'
    )

    # Map common UniProt MOD_RES /note values to (PTM type, residue)
    _MODRES_TYPE_MAP = {
        'phosphoserine': ('Phosphorylation', 'S'),
        'phosphothreonine': ('Phosphorylation', 'T'),
        'phosphotyrosine': ('Phosphorylation', 'Y'),
        'n6-acetyllysine': ('Acetylation', 'K'),
        'n-acetylalanine': ('Acetylation', 'A'),
        'n-acetylmethionine': ('Acetylation', 'M'),
        'n-acetylserine': ('Acetylation', 'S'),
        'n-acetylglycine': ('Acetylation', 'G'),
        'n-acetylthreonine': ('Acetylation', 'T'),
        'n-acetylaspartate': ('Acetylation', 'D'),
        'n-acetylglutamate': ('Acetylation', 'E'),
        'n-acetylvaline': ('Acetylation', 'V'),
        'omega-n-methylarginine': ('Methylation', 'R'),
        'asymmetric dimethylarginine': ('Methylation', 'R'),
        'symmetric dimethylarginine': ('Methylation', 'R'),
        'n6-methyllysine': ('Methylation', 'K'),
        'n6,n6-dimethyllysine': ('Methylation', 'K'),
        'n6,n6,n6-trimethyllysine': ('Methylation', 'K'),
        'sulfotyrosine': ('Sulfation', 'Y'),
        'citrulline': ('Citrullination', 'R'),
        'deamidated asparagine': ('Deamidation', 'N'),
        'deamidated glutamine': ('Deamidation', 'Q'),
        'pyrrolidone carboxylic acid': ('Pyroglutamylation', 'Q'),
    }

    def _parse_modres_note(self, note):
        """
        Parse a UniProt MOD_RES /note value to extract PTM type and residue.

        Returns:
            (ptm_type, residue) tuple
        """
        # Remove trailing qualifiers like "; alternate", "; by autocatalysis"
        note_core = note.lower().split(';')[0].strip()

        # Direct lookup
        if note_core in self._MODRES_TYPE_MAP:
            return self._MODRES_TYPE_MAP[note_core]

        # Fuzzy matching for common patterns
        if 'phospho' in note_core:
            return ('Phosphorylation', '')
        if 'acetyl' in note_core:
            return ('Acetylation', '')
        if 'methyl' in note_core:
            return ('Methylation', '')
        if 'ubiquitin' in note_core:
            return ('Ubiquitination', '')
        if 'sumo' in note_core:
            return ('SUMOylation', '')
        if 'hydroxy' in note_core:
            return ('Hydroxylation', '')
        if 'glyco' in note_core:
            return ('Glycosylation', '')
        if 'sulfo' in note_core:
            return ('Sulfation', '')
        if 'citrullin' in note_core:
            return ('Citrullination', '')
        if 'nitro' in note_core:
            return ('S-Nitrosylation', '')
        if 'palmit' in note_core:
            return ('Palmitoylation', '')
        if 'myrist' in note_core:
            return ('Myristoylation', '')

        # Unrecognised -- keep the raw note as type
        return (note.split(';')[0].strip(), '')

    def _extract_pmids_from_evidence(self, evidence):
        """Extract PubMed IDs from UniProt evidence strings."""
        if not evidence:
            return ''
        pmids = re.findall(r'PubMed:(\d+)', evidence)
        return '|'.join(sorted(set(pmids)))

    def _parse_uniprot_modres(self, header_line, fh, filename=""):
        """
        Parse UniProt MOD_RES flat annotation format.

        Expected columns: Entry, Entry Name, Gene Names, Protein names,
                          Modified residue

        The 'Modified residue' column contains concatenated MOD_RES
        annotation strings, each of the form:
            MOD_RES <position>; /note="<modification>"; /evidence="..."
        """
        header = header_line.split('\t')
        count = 0

        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                continue

            row = dict(zip(header, parts))

            acc = (row.get('Entry', '') or '').strip()
            if not acc:
                continue

            modres_col = (
                row.get('Modified residue', '')
                or row.get('Modified_residue', '')
                or ''
            ).strip()

            if not modres_col:
                continue

            # Parse all MOD_RES entries from the concatenated string
            for match in self._MODRES_PATTERN.finditer(modres_col):
                position_str = match.group(1)
                note = match.group(2)
                evidence = match.group(3) or ''

                ptm_type, residue = self._parse_modres_note(note)
                pmids = self._extract_pmids_from_evidence(evidence)

                self.ptms.append({
                    'acc': acc,
                    'position': position_str,
                    'ptm_type': ptm_type,
                    'residue': residue,
                    'pmids': pmids,
                })
                count += 1

        if count > 0:
            logger.info(
                f"dbPTM: Parsed {count} records from {filename} "
                "(UniProt MOD_RES format)"
            )

    @staticmethod
    def _normalise_ptm_type(raw):
        """Map common dbPTM type strings to canonical names."""
        if not raw:
            return ''
        lower = raw.lower().strip()
        mapping = {
            'phosphorylation': 'Phosphorylation',
            'acetylation': 'Acetylation',
            'ubiquitination': 'Ubiquitination',
            'ubiquitylation': 'Ubiquitination',
            'sumoylation': 'SUMOylation',
            'methylation': 'Methylation',
            'glycosylation': 'Glycosylation',
            'n-linked glycosylation': 'N-Glycosylation',
            'o-linked glycosylation': 'O-Glycosylation',
            's-nitrosylation': 'S-Nitrosylation',
            'succinylation': 'Succinylation',
            'malonylation': 'Malonylation',
            'palmitoylation': 'Palmitoylation',
            'myristoylation': 'Myristoylation',
            'hydroxylation': 'Hydroxylation',
            'crotonylation': 'Crotonylation',
            'neddylation': 'Neddylation',
            'sulfation': 'Sulfation',
            'sulphation': 'Sulfation',
            'citrullination': 'Citrullination',
            'nitration': 'Nitration',
            'oxidation': 'Oxidation',
            'glutathionylation': 'Glutathionylation',
            'formylation': 'Formylation',
        }
        return mapping.get(lower, raw.strip())

    def get_nodes(self):
        """
        No new nodes -- PTM sites are represented as edge properties linking
        existing protein nodes to PTM target identifiers.
        """
        logger.info("dbPTM: No dedicated nodes (uses existing protein nodes)")
        return
        yield  # noqa: unreachable -- makes this a generator

    def get_edges(self):
        """
        Generate ProteinHasPTM edges for post-translational modifications.

        Each edge connects a UniProt protein accession to a PTM site
        identifier of the form dbptm:<acc>_<position>_<ptm_type>.
        Deduplicates by (accession, position, ptm_type).

        Yields:
            (edge_id, source_id, target_id, label, properties)
        """
        logger.info("dbPTM: Generating ProteinHasPTM edges...")
        seen = set()
        count = 0

        for ptm in self.ptms:
            acc = ptm['acc']
            position = ptm['position']
            ptm_type = ptm['ptm_type']

            # Deduplicate by accession + position + PTM type
            key = (acc, position, ptm_type)
            if key in seen:
                continue
            seen.add(key)

            # Build target identifier
            type_tag = ptm_type.replace(' ', '') if ptm_type else 'PTM'
            residue = ptm['residue']
            if residue and position:
                site_id = f"dbptm:{acc}_{type_tag}_{residue}{position}"
            elif position:
                site_id = f"dbptm:{acc}_{type_tag}_{position}"
            else:
                site_id = f"dbptm:{acc}_{type_tag}_unknown"

            # Normalise PMIDs
            pmids_raw = ptm['pmids']
            if pmids_raw:
                pmids_clean = '|'.join(
                    p.strip()
                    for p in pmids_raw.replace(';', ',').split(',')
                    if p.strip()
                )
            else:
                pmids_clean = ''

            props = {
                'site': self._sanitize(position),
                'ptm_type': self._sanitize(ptm_type) if ptm_type else 'Unknown',
                'residue': self._sanitize(residue),
                'score': 0,
                'enzymes': '',
                'pmids': self._sanitize(pmids_clean),
                'source': 'dbPTM',
            }

            yield (None, acc, site_id, "ProteinHasPTM", props)
            count += 1

        logger.info(f"dbPTM: Generated {count} ProteinHasPTM edges")
