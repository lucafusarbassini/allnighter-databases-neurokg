"""
CellTalkDB Adapter for BioCypher.

Loads ligand-receptor interaction pairs for cell-cell communication from
CellTalkDB and generates:
- LigandReceptorInteraction edges (reuses the LIANA schema type)

CellTalkDB provides literature-supported ligand-receptor pairs curated from
experimental evidence. Gene symbols are used as identifiers, consistent with
the Gene node type used in the LIANA adapter.

Expected data directory: template_package/data/celltalkdb/
Expected file: celltalkdb_human.tsv, human_lr_pair.txt, or human_lr_pair.csv
"""

from pathlib import Path
from biocypher._logger import logger


class CellTalkDBAdapter:
    def __init__(self, data_dir="template_package/data/celltalkdb"):
        self.data_dir = Path(data_dir)
        self.lr_pairs = []
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load CellTalkDB ligand-receptor pair data."""
        if not self.data_dir.exists():
            logger.warning(
                f"CellTalkDB: Data directory not found: {self.data_dir}")
            return

        # Scan for all plausible data files
        candidates = (list(self.data_dir.glob("*.csv"))
                      + list(self.data_dir.glob("*.tsv"))
                      + list(self.data_dir.glob("*.txt")))

        if not candidates:
            logger.warning(
                f"CellTalkDB: No data files found in {self.data_dir}. "
                "Expected celltalkdb_human.tsv or human_lr_pair.txt"
            )
            return

        files_parsed = 0
        for fpath in candidates:
            try:
                with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
                    first_line = f.readline()
                    # Skip HTML error pages, 404 responses, or empty files
                    if first_line.strip().startswith('<') or \
                            first_line.strip().startswith('404') or \
                            len(first_line.strip()) < 10:
                        logger.warning(
                            f"CellTalkDB: {fpath.name} appears invalid "
                            "(likely download error), skipping"
                        )
                        continue
                    f.seek(0)
                    sep = ',' if fpath.suffix == '.csv' else '\t'
                    self._parse_file(f, sep)
                    files_parsed += 1
            except Exception as e:
                logger.warning(f"CellTalkDB: Error reading {fpath}: {e}")

        if files_parsed == 0:
            logger.warning("CellTalkDB: No valid data files could be parsed")

        logger.info(
            f"CellTalkDB: Loaded {len(self.lr_pairs)} ligand-receptor pairs"
        )

    def _parse_file(self, fh, sep='\t'):
        """
        Parse a CellTalkDB data file from a file handle.
        Detects column names flexibly and extracts ligand-receptor pairs.
        """
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue

            parts = line.split(sep)

            # First non-empty line is the header
            if header is None:
                header = [h.strip().strip('"') for h in parts]
                continue

            if len(parts) < 2:
                continue

            row = dict(zip(header, [p.strip().strip('"') for p in parts]))

            # Resolve ligand column (try known aliases)
            ligand = (row.get('ligand_gene_symbol')
                      or row.get('ligand_gene')
                      or row.get('ligand')
                      or row.get('source_genesymbol')
                      or row.get('source')
                      or '').strip()

            # Resolve receptor column (try known aliases)
            receptor = (row.get('receptor_gene_symbol')
                        or row.get('receptor_gene')
                        or row.get('receptor')
                        or row.get('target_genesymbol')
                        or row.get('target')
                        or '').strip()

            if not ligand or not receptor:
                continue

            # Extract optional fields
            lr_pair = self._sanitize(
                row.get('lr_pair', row.get('pair', '')))
            evidence = self._sanitize(
                row.get('evidence', row.get('lr_evidence', '')))
            species = self._sanitize(
                row.get('species', row.get('organism', '')))

            self.lr_pairs.append({
                'ligand': ligand,
                'receptor': receptor,
                'lr_pair': lr_pair,
                'evidence': evidence,
                'species': species if species else 'Human',
            })

    def get_nodes(self):
        """
        CellTalkDB does not introduce new node types. Gene nodes referenced
        by ligand and receptor symbols are expected to exist from other
        adapters (e.g. LIANA, HGNC). This method yields nothing.
        """
        logger.info("CellTalkDB: No new nodes (uses existing Gene nodes)")
        return
        yield  # makes this a generator

    def get_edges(self):
        """
        Generate LigandReceptorInteraction edges from CellTalkDB data.
        Uses gene symbols as source/target IDs consistent with the Gene
        node type. Deduplicates identical ligand-receptor pairs.
        Yields: (id, source, target, label, properties)
        """
        logger.info(
            f"CellTalkDB: Generating edges from "
            f"{len(self.lr_pairs)} interactions..."
        )
        count = 0
        seen = set()

        for pair in self.lr_pairs:
            ligand = pair['ligand']
            receptor = pair['receptor']

            # Deduplicate: same ligand-receptor pair appears only once
            pair_key = (ligand, receptor)
            if pair_key in seen:
                continue
            seen.add(pair_key)

            props = {
                'lr_pair': pair['lr_pair'],
                'evidence': pair['evidence'],
                'species': ['Homo sapiens'],
                'source': 'CellTalkDB',
            }

            yield (
                None,
                ligand,
                receptor,
                "LigandReceptorInteraction",
                props,
            )
            count += 1

        logger.info(
            f"CellTalkDB: Generated {count} LigandReceptorInteraction edges"
        )
