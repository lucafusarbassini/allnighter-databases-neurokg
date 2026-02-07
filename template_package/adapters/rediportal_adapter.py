"""
REDIportal Adapter for BioCypher.

Loads REDIportal v3 human RNA editing sites (hg38) and generates:
- RNAEditingSite nodes (A-to-I RNA editing events)

REDIportal is the largest database of RNA editing in humans,
cataloging millions of A-to-I editing sites from multiple tissues.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class REDIportalAdapter:
    def __init__(self, data_dir="template_package/data/rediportal"):
        self.data_dir = Path(data_dir)
        self.sites = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load REDIportal RNA editing sites."""
        path = self.data_dir / 'rediportal_hg38_v3.txt.gz'
        if not path.exists():
            logger.warning("REDIportal: data file not found")
            return

        logger.info("REDIportal: Loading RNA editing sites...")
        count = 0
        seen = set()
        max_sites = 500000  # Cap for memory management

        try:
            with gzip.open(path, 'rt', encoding='utf-8') as f:
                header = None
                for line in f:
                    parts = line.strip().split('\t')
                    if header is None:
                        header = parts
                        continue

                    if len(parts) < 30:
                        continue

                    accession = parts[0].strip()
                    region = parts[1].strip()      # chromosome
                    position = parts[2].strip()     # genomic position
                    ref = parts[3].strip()          # reference base
                    ed = parts[4].strip()           # edited base
                    strand = parts[5].strip()
                    db_source = parts[6].strip()    # A=ATLAS, D=DARNED
                    repeat_type = parts[7].strip()  # ALU, nonALU, nonREP
                    dbsnp = parts[8].strip()

                    # Gene annotation from GENCODE
                    func_region = parts[10].strip() if len(parts) > 10 else ''
                    gene = parts[11].strip() if len(parts) > 11 else ''

                    # Tissue/sample counts
                    n_samples = parts[23].strip() if len(parts) > 23 else '0'
                    n_tissues = parts[24].strip() if len(parts) > 24 else '0'

                    # Gene IDs
                    gene_id = parts[29].strip() if len(parts) > 29 else ''

                    if not accession or not region:
                        continue

                    # Deduplicate by accession
                    if accession in seen:
                        continue
                    seen.add(accession)

                    try:
                        n_samples_int = int(n_samples)
                    except ValueError:
                        n_samples_int = 0

                    try:
                        n_tissues_int = int(n_tissues)
                    except ValueError:
                        n_tissues_int = 0

                    # Prioritize well-supported sites (observed in multiple samples)
                    # For the capped set, only include sites seen in 2+ samples
                    if count >= 100000 and n_samples_int < 2:
                        continue

                    self.sites.append({
                        'accession': accession,
                        'chromosome': region,
                        'position': position,
                        'ref': ref,
                        'ed': ed,
                        'strand': strand,
                        'db_source': db_source,
                        'repeat_type': repeat_type,
                        'dbsnp': dbsnp,
                        'func_region': func_region,
                        'gene': gene,
                        'gene_id': gene_id,
                        'n_samples': n_samples_int,
                        'n_tissues': n_tissues_int,
                    })
                    count += 1

                    if count >= max_sites:
                        break

        except EOFError:
            logger.warning(f"REDIportal: Truncated gzip, loaded {count} sites")

        logger.info(f"REDIportal: Loaded {count} RNA editing sites")

    def get_nodes(self):
        """
        Generate RNAEditingSite nodes.
        Yields: (id, label, properties)
        """
        logger.info("REDIportal: Generating nodes...")
        count = 0

        for site in self.sites:
            props = {
                'chromosome': site['chromosome'],
                'position': site['position'],
                'ref_base': site['ref'],
                'edited_base': site['ed'],
                'strand': site['strand'],
                'repeat_type': self._sanitize(site['repeat_type']),
                'functional_region': self._sanitize(site['func_region']),
                'gene': self._sanitize(site['gene']),
                'gene_id': self._sanitize(site['gene_id']),
                'n_samples': site['n_samples'],
                'n_tissues': site['n_tissues'],
                'dbsnp': self._sanitize(site['dbsnp']),
                'source': 'REDIportal_v3',
            }

            yield (site['accession'], "RNAEditingSite", props)
            count += 1

        logger.info(f"REDIportal: Generated {count} RNAEditingSite nodes")

    def get_edges(self):
        """No edges."""
        logger.info("REDIportal: No edges to generate")
        return iter([])
