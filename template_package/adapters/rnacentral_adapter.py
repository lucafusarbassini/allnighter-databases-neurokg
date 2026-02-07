"""
RNAcentral Adapter for BioCypher.

Loads RNAcentral human ncRNA cross-reference mappings and generates:
- NcRNA nodes (non-coding RNA entities with cross-database mappings)

RNAcentral is a comprehensive database of non-coding RNA sequences
that provides a unified identifier space across ncRNA databases.
"""

from pathlib import Path
from biocypher._logger import logger


class RNAcentralAdapter:
    def __init__(self, data_dir="template_package/data/rnacentral"):
        self.data_dir = Path(data_dir)
        self.ncrnas = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load RNAcentral HGNC mappings for human ncRNAs."""
        path = self.data_dir / 'hgnc_mappings.tsv'
        if not path.exists():
            logger.warning("RNAcentral: HGNC mappings not found")
            return

        logger.info("RNAcentral: Loading HGNC ncRNA mappings...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                urs_id = parts[0].strip()
                db = parts[1].strip()
                ext_id = parts[2].strip()
                taxon = parts[3].strip()
                rna_type = parts[4].strip()
                gene_name = parts[5].strip()

                if taxon != '9606':
                    continue

                if urs_id not in self.ncrnas:
                    self.ncrnas[urs_id] = {
                        'urs_id': urs_id,
                        'rna_type': rna_type,
                        'gene_name': gene_name,
                        'hgnc_id': ext_id,
                        'databases': set(),
                    }
                self.ncrnas[urs_id]['databases'].add(db)
                count += 1

        # Also load miRBase mappings
        mirbase_path = self.data_dir / 'mirbase_mappings.tsv'
        if mirbase_path.exists():
            logger.info("RNAcentral: Loading miRBase mappings...")
            mir_count = 0
            with open(mirbase_path, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue

                    urs_id = parts[0].strip()
                    db = parts[1].strip()
                    ext_id = parts[2].strip()
                    taxon = parts[3].strip()
                    rna_type = parts[4].strip()
                    name = parts[5].strip()

                    if taxon != '9606':
                        continue

                    if urs_id not in self.ncrnas:
                        self.ncrnas[urs_id] = {
                            'urs_id': urs_id,
                            'rna_type': rna_type,
                            'gene_name': name,
                            'hgnc_id': '',
                            'databases': set(),
                        }
                    self.ncrnas[urs_id]['databases'].add(db)
                    mir_count += 1

            logger.info(f"RNAcentral: Added {mir_count} miRBase entries")

        logger.info(f"RNAcentral: Total {len(self.ncrnas)} unique human ncRNAs")

    def get_nodes(self):
        """
        Generate NcRNA nodes.
        Yields: (id, label, properties)
        """
        logger.info("RNAcentral: Generating nodes...")
        count = 0

        for urs_id, ncrna in self.ncrnas.items():
            props = {
                'gene_name': self._sanitize(ncrna['gene_name']),
                'rna_type': ncrna['rna_type'],
                'hgnc_id': ncrna['hgnc_id'],
                'databases': ';'.join(sorted(ncrna['databases'])),
                'source': 'RNAcentral',
            }

            yield (f"RNAcentral:{urs_id}", "NcRNA", props)
            count += 1

        logger.info(f"RNAcentral: Generated {count} NcRNA nodes")

    def get_edges(self):
        """No edges."""
        logger.info("RNAcentral: No edges to generate")
        return iter([])
