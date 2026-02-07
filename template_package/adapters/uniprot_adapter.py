"""
UniProt Adapter for BioCypher.

Loads UniProt Swiss-Prot reviewed human proteome and generates:
- UniProtProtein nodes (curated human protein annotations)

UniProt provides comprehensive, high-quality protein sequence and
functional information for all known human proteins.
"""

from pathlib import Path
from biocypher._logger import logger


class UniProtAdapter:
    def __init__(self, data_dir="template_package/data/uniprot"):
        self.data_dir = Path(data_dir)
        self.proteins = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load UniProt Swiss-Prot human proteome."""
        path = self.data_dir / 'human_swissprot.tsv'
        if not path.exists():
            logger.warning("UniProt: human proteome data not found")
            return

        logger.info("UniProt: Loading Swiss-Prot human proteome...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 4:
                    continue

                entry = parts[0].strip()
                gene_names = parts[1].strip() if len(parts) > 1 else ''
                protein_name = parts[2].strip() if len(parts) > 2 else ''
                length = parts[3].strip() if len(parts) > 3 else '0'
                subcellular = parts[4].strip() if len(parts) > 4 else ''
                function = parts[5].strip() if len(parts) > 5 else ''
                go_ids = parts[6].strip() if len(parts) > 6 else ''
                tissue = parts[8].strip() if len(parts) > 8 else ''

                if not entry:
                    continue

                try:
                    length_int = int(length)
                except ValueError:
                    length_int = 0

                primary_gene = gene_names.split()[0] if gene_names else ''

                self.proteins.append({
                    'entry': entry,
                    'gene_name': primary_gene,
                    'protein_name': protein_name[:200],
                    'length': length_int,
                    'subcellular': subcellular[:300],
                    'function': function[:300],
                    'go_ids': go_ids,
                    'tissue': tissue[:200],
                })
                count += 1

        logger.info(f"UniProt: Loaded {count} Swiss-Prot human proteins")

    def get_nodes(self):
        """
        Generate UniProtProtein nodes.
        Yields: (id, label, properties)
        """
        logger.info("UniProt: Generating nodes...")
        count = 0

        for prot in self.proteins:
            props = {
                'gene_name': self._sanitize(prot['gene_name']),
                'protein_name': self._sanitize(prot['protein_name']),
                'length': prot['length'],
                'subcellular_location': self._sanitize(prot['subcellular']),
                'function': self._sanitize(prot['function']),
                'tissue_specificity': self._sanitize(prot['tissue']),
                'source': 'UniProt_SwissProt',
            }

            yield (prot['entry'], "UniProtProtein", props)
            count += 1

        logger.info(f"UniProt: Generated {count} UniProtProtein nodes")

    def get_edges(self):
        """No edges."""
        logger.info("UniProt: No edges to generate")
        return iter([])
