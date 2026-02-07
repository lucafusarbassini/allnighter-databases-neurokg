"""
ClinVar Adapter for BioCypher.

Loads ClinVar variant-disease associations and generates:
- ClinVarVariant nodes (clinically significant genetic variants)

ClinVar is a public archive of reports of the relationships
among human variations and phenotypes, with supporting evidence.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class ClinVarAdapter:
    def __init__(self, data_dir="template_package/data/clinvar"):
        self.data_dir = Path(data_dir)
        self.variants = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ClinVar variant summary data."""
        path = self.data_dir / 'variant_summary.txt.gz'
        if not path.exists():
            logger.warning("ClinVar: variant summary not found")
            return

        logger.info("ClinVar: Loading variant summary (GRCh38 pathogenic/likely pathogenic)...")
        count = 0
        seen = set()

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 20:
                    continue

                allele_id = parts[0]
                var_type = parts[1]
                name = parts[2]
                gene_id = parts[3]
                gene_symbol = parts[4]
                clin_sig = parts[6]
                assembly = parts[16] if len(parts) > 16 else ''
                chromosome = parts[18] if len(parts) > 18 else ''
                start = parts[19] if len(parts) > 19 else ''
                phenotype_list = parts[13] if len(parts) > 13 else ''
                review_status = parts[24] if len(parts) > 24 else ''

                # Only GRCh38 and pathogenic/likely pathogenic
                if assembly != 'GRCh38':
                    continue
                if 'athogenic' not in clin_sig:
                    continue

                # Deduplicate
                if allele_id in seen:
                    continue
                seen.add(allele_id)

                try:
                    start_int = int(start)
                except ValueError:
                    start_int = 0

                self.variants.append({
                    'allele_id': allele_id,
                    'type': var_type,
                    'name': name[:200],
                    'gene_symbol': gene_symbol,
                    'clinical_significance': clin_sig,
                    'chromosome': chromosome,
                    'start': start_int,
                    'phenotype': phenotype_list[:200],
                    'review_status': review_status,
                })
                count += 1

                if count >= 500000:
                    break

        logger.info(f"ClinVar: Loaded {count} pathogenic/likely pathogenic variants")

    def get_nodes(self):
        """
        Generate ClinVarVariant nodes.
        Yields: (id, label, properties)
        """
        logger.info("ClinVar: Generating nodes...")
        count = 0

        for var in self.variants:
            props = {
                'type': var['type'],
                'name': self._sanitize(var['name']),
                'gene_symbol': self._sanitize(var['gene_symbol']),
                'clinical_significance': self._sanitize(var['clinical_significance']),
                'chromosome': var['chromosome'],
                'start': var['start'],
                'phenotype': self._sanitize(var['phenotype']),
                'review_status': self._sanitize(var['review_status']),
                'source': 'ClinVar',
            }

            yield (f"ClinVar:{var['allele_id']}", "ClinVarVariant", props)
            count += 1

        logger.info(f"ClinVar: Generated {count} ClinVarVariant nodes")

    def get_edges(self):
        """No edges."""
        logger.info("ClinVar: No edges to generate")
        return iter([])
