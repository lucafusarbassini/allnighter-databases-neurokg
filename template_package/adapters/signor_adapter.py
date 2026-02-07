"""
SIGNOR (Signaling Network Open Resource) Adapter for BioCypher.

Loads SIGNOR signaling network interactions and generates:
- SignalingInteraction edges (protein -> protein signaling)

SIGNOR is a database of causal relationships between biological
entities, manually annotated from the literature, focusing on
human signaling pathways.
"""

import csv
from pathlib import Path
from biocypher._logger import logger

# SIGNOR TSV column names.
# The bulk-download format has 28 columns: the standard 25 documented
# columns plus two additional empty/internal columns between
# MODIFICATIONB and PMID, and a SCORE column at the end.
SIGNOR_COLUMNS = [
    'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA',
    'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
    'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
    'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
    'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
    'MODIFICATIONA', 'MODIFICATIONB',
    '_RESERVED1', '_RESERVED2',
    'PMID', 'DIRECT', 'NOTES', 'ANNOTATOR',
    'SENTENCE', 'SIGNOR_ID', 'SCORE',
]


class SIGNORAdapter:
    def __init__(self, data_dir="template_package/data/signor"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _is_html(self, filepath):
        """Check if a file is HTML rather than TSV data."""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                first_line = f.readline(512).strip()
                return first_line.startswith('<!') or first_line.lower().startswith('<html')
        except Exception:
            return False

    def _find_data_file(self):
        """Find the SIGNOR TSV data file, trying known filenames."""
        candidates = [
            'signor_all_data_23_11_04.tsv',
            'signor_human.tsv',
            'signor_all_data.tsv',
        ]

        # First try known filenames
        for name in candidates:
            path = self.data_dir / name
            if path.exists() and not self._is_html(path):
                return path

        # Fall back to any TSV file in the directory that is not HTML
        if self.data_dir.exists():
            for entry in sorted(self.data_dir.iterdir()):
                if entry.suffix == '.tsv' and entry.is_file():
                    if not self._is_html(entry):
                        return entry

        return None

    def _has_header(self, filepath):
        """Check whether a TSV file has a SIGNOR header row."""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                fields = first_line.split('\t')
                return 'ENTITYA' in fields or 'IDA' in fields
        except Exception:
            return False

    def _load_data(self):
        """Load SIGNOR signaling interactions from TSV."""
        if not self.data_dir.exists():
            logger.warning("SIGNOR: data directory not found")
            return

        path = self._find_data_file()
        if path is None:
            logger.warning("SIGNOR: no valid TSV data file found (files may be HTML)")
            return

        logger.info(f"SIGNOR: Loading interactions from {path.name}...")
        count = 0
        seen = set()

        has_header = self._has_header(path)

        with open(path, 'r', encoding='utf-8') as f:
            if has_header:
                reader = csv.DictReader(f, delimiter='\t')
            else:
                # Supply standard SIGNOR column names for headerless files;
                # extra trailing columns beyond the 25 standard ones are
                # collected under the restkey and ignored.
                reader = csv.DictReader(
                    f, fieldnames=SIGNOR_COLUMNS,
                    delimiter='\t', restkey='_extra',
                )

            for row in reader:
                entity_a = (row.get('ENTITYA') or '').strip()
                entity_b = (row.get('ENTITYB') or '').strip()
                id_a = (row.get('IDA') or '').strip()
                id_b = (row.get('IDB') or '').strip()
                type_a = (row.get('TYPEA') or '').strip()
                type_b = (row.get('TYPEB') or '').strip()
                db_a = (row.get('DATABASEA') or '').strip()
                db_b = (row.get('DATABASEB') or '').strip()
                effect = (row.get('EFFECT') or '').strip()
                mechanism = (row.get('MECHANISM') or '').strip()
                residue = (row.get('RESIDUE') or '').strip()
                tax_id = (row.get('TAX_ID') or '').strip()
                direct = (row.get('DIRECT') or '').strip()
                pmid = (row.get('PMID') or '').strip()
                signor_id = (row.get('SIGNOR_ID') or '').strip()
                cell_data = (row.get('CELL_DATA') or '').strip()
                tissue_data = (row.get('TISSUE_DATA') or '').strip()
                sentence = (row.get('SENTENCE') or '').strip()

                # Filter for human interactions (TAX_ID 9606)
                if tax_id and tax_id != '9606':
                    continue

                if not id_a or not id_b:
                    continue

                # Build source and target IDs, preferring UniProt identifiers
                source_id = id_a
                if db_a.upper() == 'UNIPROT':
                    source_id = f"uniprot:{id_a}"
                elif db_a:
                    source_id = f"{db_a}:{id_a}"

                target_id = id_b
                if db_b.upper() == 'UNIPROT':
                    target_id = f"uniprot:{id_b}"
                elif db_b:
                    target_id = f"{db_b}:{id_b}"

                # Deduplicate by source-target-mechanism
                key = (source_id, target_id, mechanism)
                if key in seen:
                    continue
                seen.add(key)

                self.interactions.append({
                    'source_id': source_id,
                    'target_id': target_id,
                    'entity_a': entity_a,
                    'entity_b': entity_b,
                    'type_a': type_a,
                    'type_b': type_b,
                    'effect': effect,
                    'mechanism': mechanism,
                    'residue': residue,
                    'direct': direct,
                    'pmid': pmid,
                    'signor_id': signor_id,
                    'cell_data': cell_data,
                    'tissue_data': tissue_data,
                    'sentence': sentence,
                })
                count += 1

        logger.info(f"SIGNOR: Loaded {count} human signaling interactions")

    def get_nodes(self):
        """No new nodes -- proteins are expected from UniProt adapter."""
        logger.info("SIGNOR: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate SignalingInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("SIGNOR: Generating edges...")
        count = 0

        for ix in self.interactions:
            props = {
                'entity_a': self._sanitize(ix['entity_a']),
                'entity_b': self._sanitize(ix['entity_b']),
                'type_a': self._sanitize(ix['type_a']),
                'type_b': self._sanitize(ix['type_b']),
                'effect': self._sanitize(ix['effect']),
                'mechanism': self._sanitize(ix['mechanism']),
                'residue': self._sanitize(ix['residue']),
                'is_direct': ix['direct'],
                'pmid': self._sanitize(ix['pmid']),
                'signor_id': self._sanitize(ix['signor_id']),
                'cell_data': self._sanitize(ix['cell_data']),
                'tissue_data': self._sanitize(ix['tissue_data']),
                'sentence': self._sanitize(
                    ix['sentence'][:300] if ix['sentence'] else ''
                ),
                'source': 'SIGNOR',
            }

            yield (
                None,
                ix['source_id'],
                ix['target_id'],
                "SignalingInteraction",
                props
            )
            count += 1

        logger.info(f"SIGNOR: Generated {count} SignalingInteraction edges")
