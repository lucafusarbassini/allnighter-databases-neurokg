"""
sORFs.org (Small Open Reading Frames) Adapter for BioCypher.

Loads small ORF data from ribosome profiling.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class SORFsAdapter:
    def __init__(self, data_dir="template_package/data/sorfs"):
        self.data_dir = Path(data_dir)
        self.sorfs = []
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
            logger.warning("sORFs: data directory not found")
            return
        candidates = (list(self.data_dir.glob("*.tsv"))
                      + list(self.data_dir.glob("*.txt"))
                      + list(self.data_dir.glob("*.csv"))
                      + list(self.data_dir.glob("*.tsv.gz")))
        for fpath in candidates:
            try:
                if str(fpath).endswith('.gz'):
                    with gzip.open(fpath, 'rt', errors='replace') as f:
                        first = f.readline()
                        if first.startswith('<'):
                            continue
                        f.seek(0)
                        self._parse_file(f)
                else:
                    with open(fpath, 'r', errors='replace') as f:
                        first = f.readline()
                        if first.startswith('<'):
                            continue
                        f.seek(0)
                        self._parse_file(f)
            except Exception as e:
                logger.warning(f"sORFs: Error reading {fpath}: {e}")
        logger.info(f"sORFs: Loaded {len(self.sorfs)} small ORFs")

    def _parse_file(self, fh):
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            if len(parts) < 3:
                continue
            row = dict(zip(header, parts))
            self.sorfs.append(row)

    def get_nodes(self):
        logger.info("sORFs: Generating AlternativeProtein nodes...")
        count = 0
        for sorf in self.sorfs:
            sorf_id = sorf.get('sORF_ID', sorf.get('id', f"sorf_{count}"))
            gene = sorf.get('gene', sorf.get('gene_symbol', ''))
            props = {
                'protein_type': 'sORF',
                'gene_name': self._sanitize(gene),
                'sequence_length': int(sorf.get('length', 0)) if sorf.get('length', '').isdigit() else 0,
                'transcript_accessions': self._sanitize(sorf.get('transcript', '')),
                'source': 'sORFs.org',
            }
            yield (f"sorf:{sorf_id}", "AlternativeProtein", props)
            count += 1
        logger.info(f"sORFs: Generated {count} AlternativeProtein nodes")

    def get_edges(self):
        logger.info("sORFs: No edges")
        return
        yield
