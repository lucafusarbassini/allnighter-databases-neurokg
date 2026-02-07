"""
Phospho.ELM Adapter for BioCypher.

Loads Phospho.ELM phosphorylation site data and generates:
- ProteinHasPTM edges (protein -> phosphorylation site with kinase and evidence)

Phospho.ELM contains experimentally verified phosphorylation sites (Ser/Thr/Tyr)
in eukaryotic proteins, manually curated from the scientific literature.

Expected data directory: template_package/data/phosphoelm/
Expected file: phosphoELM_vertebrate_latest.dump or phosphoELM_all_latest.dump
  (tab-separated with columns: acc, sequence, position, code, pmids, kinases,
   source, species, entry_date)

Reference:
  Diella et al., Phospho.ELM: a database of phosphorylation sites--update 2008.
  Nucleic Acids Res. 2008 Jan;36(Database issue):D240-4.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class PhosphoELMAdapter:
    def __init__(self, data_dir="template_package/data/phosphoelm"):
        self.data_dir = Path(data_dir)
        self.sites = []
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
        Locate and load PhosphoELM dump files.

        Searches the data directory for .dump, .dump.gz, .tsv, and .txt files.
        Each candidate file is checked: if it starts with '<' (HTML), it is
        skipped, which handles the common case where the download gate returns
        an HTML form instead of the actual data.
        """
        if not self.data_dir.exists():
            logger.warning("PhosphoELM: data directory not found")
            return

        candidates = (
            list(self.data_dir.glob("*.dump"))
            + list(self.data_dir.glob("*.dump.gz"))
            + list(self.data_dir.glob("*.dump.tgz"))
            + list(self.data_dir.glob("*.tsv"))
            + list(self.data_dir.glob("*.tsv.gz"))
            + list(self.data_dir.glob("*.txt"))
        )

        if not candidates:
            logger.warning("PhosphoELM: no data files found in data directory")
            return

        for fpath in candidates:
            try:
                if str(fpath).endswith('.tgz') or str(fpath).endswith('.gz'):
                    with gzip.open(fpath, 'rt', errors='replace') as f:
                        first = f.readline()
                        if first.lstrip().startswith('<'):
                            logger.warning(
                                f"PhosphoELM: {fpath.name} appears to be "
                                "HTML (download gate?), skipping"
                            )
                            continue
                        f.seek(0)
                        self._parse_dump(f, fpath.name)
                else:
                    with open(fpath, 'r', errors='replace') as f:
                        first = f.readline()
                        if first.lstrip().startswith('<'):
                            logger.warning(
                                f"PhosphoELM: {fpath.name} appears to be "
                                "HTML (download gate?), skipping"
                            )
                            continue
                        f.seek(0)
                        self._parse_dump(f, fpath.name)
            except (gzip.BadGzipFile, OSError) as e:
                logger.warning(
                    f"PhosphoELM: Cannot read {fpath.name} (corrupt?): {e}"
                )
            except Exception as e:
                logger.warning(f"PhosphoELM: Error reading {fpath.name}: {e}")

        logger.info(
            f"PhosphoELM: Loaded {len(self.sites)} human phosphorylation sites"
        )

    def _parse_dump(self, fh, filename=""):
        """
        Parse a PhosphoELM tab-separated dump file.

        Expected columns: acc, sequence, position, code, pmids, kinases,
                          source, species, entry_date

        Filters for human entries (species contains 'Homo sapiens').
        Skips rows with missing accession or insufficient columns.
        """
        header = None
        count_before = len(self.sites)

        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Safety check for HTML content mid-file
            if line.startswith('<!DOCTYPE') or line.startswith('<html'):
                logger.warning(
                    f"PhosphoELM: {filename} contains HTML content, "
                    "stopping parse"
                )
                return

            parts = line.split('\t')

            if header is None:
                header = parts
                continue

            if len(parts) < 3:
                continue

            row = dict(zip(header, parts))

            # Extract accession -- PhosphoELM uses 'acc' as primary column
            acc = (
                row.get('acc', '')
                or row.get('ACC', '')
                or row.get('substrate', '')
                or row.get('accession', '')
            ).strip()

            if not acc:
                continue

            # Filter for human entries
            species = (
                row.get('species', '')
                or row.get('Species', '')
            ).strip()
            if species and 'Homo sapiens' not in species:
                continue

            position = (
                row.get('position', '')
                or row.get('Position', '')
            ).strip()

            # The 'code' field in PhosphoELM is the phosphorylated residue
            # (S, T, or Y), not to be confused with position
            residue = (
                row.get('code', '')
                or row.get('Code', '')
                or row.get('residue', '')
            ).strip()

            kinases = (
                row.get('kinases', '')
                or row.get('Kinases', '')
                or row.get('kinase', '')
            ).strip()

            pmids = (
                row.get('pmids', '')
                or row.get('PMIDs', '')
                or row.get('pmid', '')
            ).strip()

            source_type = (
                row.get('source', '')
                or row.get('Source', '')
            ).strip()

            sequence = (
                row.get('sequence', '')
                or row.get('Sequence', '')
            ).strip()

            entry_date = (
                row.get('entry_date', '')
                or row.get('Entry_date', '')
            ).strip()

            self.sites.append({
                'acc': acc,
                'position': position,
                'residue': residue,
                'kinases': kinases,
                'pmids': pmids,
                'source_type': source_type,
                'sequence': sequence[:15] if sequence else '',
                'entry_date': entry_date,
            })

        added = len(self.sites) - count_before
        if added > 0:
            logger.info(
                f"PhosphoELM: Parsed {added} entries from {filename}"
            )

    def get_nodes(self):
        """
        No new nodes -- phosphorylation sites are represented as edge
        properties linking existing protein nodes to PTM target identifiers.
        """
        logger.info("PhosphoELM: No dedicated nodes (uses existing protein nodes)")
        return
        yield  # noqa: unreachable -- makes this a generator

    def get_edges(self):
        """
        Generate ProteinHasPTM edges for phosphorylation sites.

        Each edge connects a UniProt protein accession to a phosphosite
        identifier of the form phosphoelm:<acc>_<residue><position>.

        Yields:
            (edge_id, source_id, target_id, label, properties)
        """
        logger.info("PhosphoELM: Generating ProteinHasPTM edges...")
        seen = set()
        count = 0

        for site in self.sites:
            acc = site['acc']
            position = site['position']
            residue = site['residue']

            # Deduplicate by accession + position + residue
            key = (acc, position, residue)
            if key in seen:
                continue
            seen.add(key)

            # Build phosphosite target identifier
            if residue and position:
                site_id = f"phosphoelm:{acc}_{residue}{position}"
            elif position:
                site_id = f"phosphoelm:{acc}_{position}"
            else:
                site_id = f"phosphoelm:{acc}_unknown"

            # Normalise kinase list (may be semicolon- or comma-separated)
            kinases_raw = site['kinases']
            if kinases_raw:
                kinases_clean = '|'.join(
                    k.strip()
                    for k in kinases_raw.replace(';', ',').split(',')
                    if k.strip()
                )
            else:
                kinases_clean = ''

            # Normalise PMIDs
            pmids_raw = site['pmids']
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
                'ptm_type': 'Phosphorylation',
                'residue': self._sanitize(residue),
                'score': 0,
                'enzymes': self._sanitize(kinases_clean),
                'pmids': self._sanitize(pmids_clean),
                'source_type': self._sanitize(site['source_type']),
                'sequence_window': self._sanitize(site['sequence']),
                'entry_date': self._sanitize(site['entry_date']),
                'source': 'PhosphoELM',
            }

            yield (None, acc, site_id, "ProteinHasPTM", props)
            count += 1

        logger.info(f"PhosphoELM: Generated {count} ProteinHasPTM edges")
