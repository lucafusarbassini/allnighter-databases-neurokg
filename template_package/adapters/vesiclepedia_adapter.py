"""
Vesiclepedia (Extracellular Vesicle Database) Adapter for BioCypher.

Loads protein and RNA content of extracellular vesicles.
- GeneInExosome edges (reuses ExoCarta schema)
"""

from pathlib import Path
from biocypher._logger import logger


class VesiclepediaAdapter:
    def __init__(self, data_dir="template_package/data/vesiclepedia"):
        self.data_dir = Path(data_dir)
        self.entries = []
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
            logger.warning("Vesiclepedia: data directory not found")
            return
        candidates = (list(self.data_dir.glob("*.tsv"))
                      + list(self.data_dir.glob("*.txt"))
                      + list(self.data_dir.glob("*.csv")))
        for fpath in candidates:
            try:
                with open(fpath, 'r', errors='replace') as f:
                    first = f.readline()
                    if first.startswith('<'):
                        continue
                    f.seek(0)
                    self._parse_file(f)
            except Exception as e:
                logger.warning(f"Vesiclepedia: Error reading {fpath}: {e}")
        logger.info(f"Vesiclepedia: Loaded {len(self.entries)} entries")

    def _get_col(self, row, *candidates):
        """Case-insensitive column lookup."""
        for c in candidates:
            if c in row:
                return row[c]
        cl = {k.lower(): v for k, v in row.items()}
        for c in candidates:
            if c.lower() in cl:
                return cl[c.lower()]
        return ''

    def _parse_file(self, fh):
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = [h.strip() for h in parts]
                continue
            if len(parts) < 3:
                continue
            row = dict(zip(header, parts))
            gene = self._get_col(row, 'GENE SYMBOL', 'Gene Symbol',
                                 'gene_symbol', 'Gene', 'Protein Name')
            species = self._get_col(row, 'SPECIES', 'Species', 'Organism')
            if species and 'Homo sapiens' not in species and 'Human' not in species:
                if 'sapiens' not in species.lower() and '9606' not in str(species):
                    continue
            if not gene:
                continue
            content_type = self._get_col(row, 'CONTENT TYPE', 'Content Type',
                                         'content_type') or 'protein'
            methods = self._get_col(row, 'METHODS', 'Vesicle Type', 'vesicle_type')
            self.entries.append({
                'gene': gene.strip(),
                'species': species,
                'content_type': content_type,
                'vesicle_type': methods,
                'tissue': self._get_col(row, 'Tissue', 'Sample'),
            })

    def get_nodes(self):
        logger.info("Vesiclepedia: No dedicated nodes")
        return
        yield

    def get_edges(self):
        logger.info("Vesiclepedia: Generating GeneInExosome edges...")
        count = 0
        seen = set()
        for entry in self.entries:
            key = (entry['gene'], entry.get('content_type', ''))
            if key in seen:
                continue
            seen.add(key)
            props = {
                'species': 'Homo sapiens',
                'content_types': self._sanitize(entry.get('content_type', '')),
                'detection_methods': self._sanitize(entry.get('vesicle_type', '')),
                'num_experiments': 1,
                'source': 'Vesiclepedia',
            }
            yield (None, entry['gene'], f"vesiclepedia:{entry['gene']}",
                   "GeneInExosome", props)
            count += 1
        logger.info(f"Vesiclepedia: Generated {count} GeneInExosome edges")
