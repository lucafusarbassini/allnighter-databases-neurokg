"""
GENCODE Adapter for BioCypher.

Loads GENCODE v46 human gene annotations from GTF and generates:
- GencodeGene nodes (comprehensive gene-level annotations)

GENCODE provides reference gene annotation for the human genome,
including protein-coding genes, lncRNAs, pseudogenes, and more.
"""

import gzip
import re
from pathlib import Path
from biocypher._logger import logger


class GENCODEAdapter:
    def __init__(self, data_dir="template_package/data/gencode"):
        self.data_dir = Path(data_dir)
        self.genes = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _parse_attributes(self, attr_str):
        """Parse GTF attribute string into a dict."""
        attrs = {}
        for match in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
            key, val = match.group(1), match.group(2)
            attrs[key] = val
        return attrs

    def _load_data(self):
        """Load GENCODE gene annotations from GTF."""
        path = self.data_dir / 'gencode.v46.annotation.gtf.gz'
        if not path.exists():
            logger.warning("GENCODE: GTF file not found")
            return

        logger.info("GENCODE: Loading gene annotations...")
        count = 0

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'gene':
                    continue

                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = self._parse_attributes(parts[8])

                gene_id = attrs.get('gene_id', '')
                gene_name = attrs.get('gene_name', '')
                gene_type = attrs.get('gene_type', '')
                level = attrs.get('level', '')
                hgnc_id = attrs.get('hgnc_id', '')

                if not gene_id:
                    continue

                self.genes.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'level': level,
                    'hgnc_id': hgnc_id,
                })
                count += 1

        logger.info(f"GENCODE: Loaded {count} gene annotations")

    def get_nodes(self):
        """
        Generate GencodeGene nodes.
        Yields: (id, label, properties)
        """
        logger.info("GENCODE: Generating nodes...")
        count = 0

        for gene in self.genes:
            props = {
                'gene_name': self._sanitize(gene['gene_name']),
                'gene_type': gene['gene_type'],
                'chromosome': gene['chromosome'],
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'annotation_level': gene['level'],
                'hgnc_id': gene['hgnc_id'],
                'source': 'GENCODE_v46',
            }

            yield (gene['gene_id'], "GencodeGene", props)
            count += 1

        logger.info(f"GENCODE: Generated {count} GencodeGene nodes")

    def get_edges(self):
        """No edges."""
        logger.info("GENCODE: No edges to generate")
        return iter([])
