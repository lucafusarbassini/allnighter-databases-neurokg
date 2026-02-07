"""
DisGeNET Adapter for BioCypher.

Loads DisGeNET gene-disease association data and generates:
- Disease nodes (human diseases with ontology IDs)
- GeneDiseaseAssociation edges (gene -> disease associations)

DisGeNET is a discovery platform containing collections of genes
and variants associated to human diseases, integrating data from
expert-curated repositories and text-mining derived associations.
"""

import csv
import gzip
from pathlib import Path
from biocypher._logger import logger


class DisGeNETAdapter:
    def __init__(self, data_dir="template_package/data/disgenet"):
        self.data_dir = Path(data_dir)
        self.diseases = {}
        self.associations = []
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
            with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
                first_line = f.readline(512).strip()
                return (
                    first_line.startswith('<!')
                    or first_line.lower().startswith('<html')
                )
        except Exception:
            return False

    def _is_gzipped(self, filepath):
        """Check if a file is actually gzip-compressed by reading magic bytes."""
        try:
            with open(filepath, 'rb') as f:
                magic = f.read(2)
                return magic == b'\x1f\x8b'
        except Exception:
            return False

    def _find_data_file(self):
        """
        Find the DisGeNET TSV data file, trying known filenames.

        Returns a tuple of (path, is_gzipped) or (None, False) if no
        valid file is found. Files that are actually HTML (common when
        the download requires authentication) are skipped.
        """
        candidates = [
            'curated_gene_disease.tsv',
            'curated_gene_disease_associations.tsv',
            'curated_gene_disease.tsv.gz',
            'curated_gene_disease_associations.tsv.gz',
            'all_gene_disease_associations.tsv',
            'all_gene_disease_associations.tsv.gz',
        ]

        for name in candidates:
            path = self.data_dir / name
            if not path.exists():
                continue

            if name.endswith('.gz'):
                if self._is_gzipped(path):
                    return path, True
                # File has .gz extension but is not actually compressed
                if not self._is_html(path):
                    return path, False
            else:
                if not self._is_html(path):
                    return path, False

        # Fall back: search for any TSV file in directory
        if self.data_dir.exists():
            for entry in sorted(self.data_dir.iterdir()):
                if not entry.is_file() or '.tsv' not in entry.name:
                    continue
                if entry.name.endswith('.gz') and self._is_gzipped(entry):
                    return entry, True
                if not self._is_html(entry):
                    return entry, False

        return None, False

    def _open_file(self, path, is_gzipped):
        """Open a file, handling gzip transparently."""
        if is_gzipped:
            return gzip.open(path, 'rt', encoding='utf-8')
        return open(path, 'r', encoding='utf-8', errors='replace')

    def _load_data(self):
        """Load DisGeNET curated gene-disease associations."""
        if not self.data_dir.exists():
            logger.warning("DisGeNET: data directory not found")
            return

        path, is_gzipped = self._find_data_file()
        if path is None:
            logger.warning(
                "DisGeNET: no valid TSV data file found (files may be HTML)"
            )
            return

        logger.info(
            f"DisGeNET: Loading gene-disease associations from {path.name}..."
        )
        count = 0
        seen = set()

        with self._open_file(path, is_gzipped) as f:
            reader = csv.DictReader(f, delimiter='\t')

            for row in reader:
                gene_id = (row.get('geneId') or '').strip()
                gene_symbol = (row.get('geneSymbol') or '').strip()
                disease_id = (row.get('diseaseId') or '').strip()
                disease_name = (row.get('diseaseName') or '').strip()
                disease_type = (row.get('diseaseType') or '').strip()
                disease_class = (row.get('diseaseClass') or '').strip()
                disease_semantic_type = (
                    row.get('diseaseSemanticType') or ''
                ).strip()
                score_str = (row.get('score') or '').strip()
                ei_str = (row.get('EI') or '').strip()
                year_initial = (row.get('YearInitial') or '').strip()
                year_final = (row.get('YearFinal') or '').strip()
                n_pmids = (row.get('NofPmids') or '').strip()
                n_snps = (row.get('NofSnps') or '').strip()
                dsi_str = (row.get('DSI') or '').strip()
                dpi_str = (row.get('DPI') or '').strip()
                source = (row.get('source') or '').strip()

                if not gene_id or not disease_id:
                    continue

                # Parse numeric fields
                try:
                    score = float(score_str)
                except (ValueError, TypeError):
                    score = 0.0

                try:
                    ei = float(ei_str)
                except (ValueError, TypeError):
                    ei = 0.0

                try:
                    pmid_count = int(n_pmids)
                except (ValueError, TypeError):
                    pmid_count = 0

                try:
                    snp_count = int(n_snps)
                except (ValueError, TypeError):
                    snp_count = 0

                # Deduplicate by gene-disease pair
                key = (gene_id, disease_id)
                if key in seen:
                    continue
                seen.add(key)

                # Register disease node
                if disease_id not in self.diseases:
                    self.diseases[disease_id] = {
                        'name': disease_name,
                        'type': disease_type,
                        'disease_class': disease_class,
                        'semantic_type': disease_semantic_type,
                    }

                self.associations.append({
                    'gene_id': gene_id,
                    'gene_symbol': gene_symbol,
                    'disease_id': disease_id,
                    'score': score,
                    'ei': ei,
                    'year_initial': year_initial,
                    'year_final': year_final,
                    'pmid_count': pmid_count,
                    'snp_count': snp_count,
                    'dsi': dsi_str,
                    'dpi': dpi_str,
                    'source': source,
                })
                count += 1

        logger.info(
            f"DisGeNET: Loaded {count} gene-disease associations, "
            f"{len(self.diseases)} diseases"
        )

    def get_nodes(self):
        """
        Generate Disease nodes.
        Yields: (id, label, properties)
        """
        logger.info("DisGeNET: Generating Disease nodes...")
        count = 0

        for disease_id, data in self.diseases.items():
            props = {
                'name': self._sanitize(data['name']),
                'disease_type': self._sanitize(data['type']),
                'disease_class': self._sanitize(
                    data['disease_class'][:300]
                    if data['disease_class']
                    else ''
                ),
                'semantic_type': self._sanitize(data['semantic_type']),
                'source': 'DisGeNET',
            }

            yield (disease_id, "Disease", props)
            count += 1

        logger.info(f"DisGeNET: Generated {count} Disease nodes")

    def get_edges(self):
        """
        Generate GeneDiseaseAssociation edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("DisGeNET: Generating edges...")
        count = 0

        for assoc in self.associations:
            # Source is NCBI gene ID, target is disease ID (UMLS CUI)
            source_id = f"ncbigene:{assoc['gene_id']}"
            target_id = assoc['disease_id']

            props = {
                'gene_symbol': self._sanitize(assoc['gene_symbol']),
                'score': assoc['score'],
                'ei': assoc['ei'],
                'year_initial': self._sanitize(assoc['year_initial']),
                'year_final': self._sanitize(assoc['year_final']),
                'pmid_count': assoc['pmid_count'],
                'snp_count': assoc['snp_count'],
                'dsi': self._sanitize(assoc['dsi']),
                'dpi': self._sanitize(assoc['dpi']),
                'disgenet_source': self._sanitize(assoc['source']),
                'source': 'DisGeNET',
            }

            yield (
                None,
                source_id,
                target_id,
                "GeneDiseaseAssociation",
                props
            )
            count += 1

        logger.info(f"DisGeNET: Generated {count} GeneDiseaseAssociation edges")
