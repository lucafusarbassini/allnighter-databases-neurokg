"""
TREP-db (Transposable Element Database) Adapter for BioCypher.

Loads repetitive/transposable element data.
"""

from pathlib import Path
from biocypher._logger import logger


class TREPdbAdapter:
    def __init__(self, data_dir="template_package/data/trepdb"):
        self.data_dir = Path(data_dir)
        self.elements = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("TREP-db: data directory not found")
            return
        for fpath in list(self.data_dir.glob("*.fasta")) + list(self.data_dir.glob("*.fa")) + list(self.data_dir.glob("*.fasta.gz")):
            try:
                self._parse_fasta(fpath)
            except Exception as e:
                logger.warning(f"TREP-db: Error reading {fpath}: {e}")
        logger.info(f"TREP-db: Loaded {len(self.elements)} elements")

    def _parse_fasta(self, fpath):
        import gzip
        opener = gzip.open if str(fpath).endswith('.gz') else open
        with opener(fpath, 'rt', errors='replace') as f:
            for line in f:
                if line.startswith('>'):
                    parts = line[1:].strip().split()
                    te_id = parts[0] if parts else ''
                    desc = ' '.join(parts[1:]) if len(parts) > 1 else ''
                    self.elements.append({
                        'id': te_id,
                        'description': desc,
                    })

    def get_nodes(self):
        logger.info("TREP-db: Generating TransposableElementFamily nodes...")
        count = 0
        for elem in self.elements:
            props = {
                'name': self._sanitize(elem['id']),
                'title': self._sanitize(elem['description']),
                'consensus_length': 0,
                'repeat_type': '',
                'repeat_subtype': '',
                'classification': '',
                'source': 'TREP-db',
            }
            yield (f"trepdb:{elem['id']}", "TransposableElementFamily", props)
            count += 1
        logger.info(f"TREP-db: Generated {count} TransposableElementFamily nodes")

    def get_edges(self):
        logger.info("TREP-db: No edges")
        return
        yield
