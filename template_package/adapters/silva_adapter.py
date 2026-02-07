"""
SILVA Adapter for BioCypher.

Loads SILVA ribosomal RNA taxonomy data. SILVA provides comprehensive,
quality-checked and regularly updated databases of aligned ribosomal
RNA (rRNA) gene sequences from bacteria, archaea, and eukaryota.

While primarily a microbiology/taxonomy resource, this is useful for
understanding the broader biological context of organisms referenced
in other databases.

Generates:
- Taxonomy nodes (taxonomic classification hierarchy)
- Taxonomy parent-child edges
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class SILVAAdapter:
    def __init__(self, data_dir="template_package/data/silva"):
        self.data_dir = Path(data_dir)
        self.taxa = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load SILVA taxonomy data."""
        # Load SSU taxonomy (most comprehensive)
        tax_path = self.data_dir / 'tax_slv_ssu_138.2.txt.gz'
        if not tax_path.exists():
            # Try alternate version
            for alt in self.data_dir.glob('tax_slv_ssu_*.txt.gz'):
                tax_path = alt
                break

        if not tax_path.exists():
            logger.warning("SILVA: No taxonomy file found")
            return

        logger.info(f"SILVA: Loading taxonomy from {tax_path.name}...")
        count = 0

        try:
            with gzip.open(tax_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    # Format: path\ttaxid\trank\tremark\trelease
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue

                    tax_path_str = parts[0].strip()
                    tax_id = parts[1].strip()
                    rank = parts[2].strip()

                    if not tax_path_str:
                        continue

                    # Extract the last taxon name from the path
                    path_parts = [p.strip() for p in tax_path_str.rstrip(';').split(';') if p.strip()]
                    if not path_parts:
                        continue

                    name = path_parts[-1]
                    full_path = ';'.join(path_parts)

                    # Determine parent
                    parent_path = ';'.join(path_parts[:-1]) if len(path_parts) > 1 else ''

                    self.taxa[full_path] = {
                        'tax_id': tax_id,
                        'name': name,
                        'rank': rank,
                        'full_path': full_path,
                        'parent_path': parent_path,
                    }
                    count += 1

        except Exception as e:
            logger.warning(f"SILVA: Error reading taxonomy: {e}")

        logger.info(f"SILVA: Loaded {count} taxonomic entries")

    def get_nodes(self):
        """
        Generate taxonomy nodes.
        Yields: (id, label, properties)
        """
        logger.info("SILVA: Generating taxonomy nodes...")
        count = 0

        for path, taxon in self.taxa.items():
            props = {
                'name': self._sanitize(taxon['name']),
                'rank': taxon['rank'],
                'full_path': self._sanitize(taxon['full_path']),
                'source': 'SILVA',
            }

            yield (f"SILVA:{taxon['tax_id']}", "TaxonomicEntity", props)
            count += 1

        logger.info(f"SILVA: Generated {count} taxonomy nodes")

    def get_edges(self):
        """
        Generate taxonomy parent-child edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("SILVA: Generating taxonomy hierarchy edges...")
        count = 0

        for path, taxon in self.taxa.items():
            parent_path = taxon['parent_path']
            if parent_path and parent_path in self.taxa:
                parent = self.taxa[parent_path]
                yield (
                    None,
                    f"SILVA:{taxon['tax_id']}",
                    f"SILVA:{parent['tax_id']}",
                    "TaxonChildOf",
                    {'source': 'SILVA'}
                )
                count += 1

        logger.info(f"SILVA: Generated {count} taxonomy hierarchy edges")
