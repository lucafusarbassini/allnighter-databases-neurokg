"""
MatrisomeDB (Extracellular Matrix) Adapter for BioCypher.

Loads MatrisomeDB proteomics data and generates:
- MatrisomeProtein nodes (ECM proteins with tissue annotations)
- ProteinDetectedInTissue edges (protein → tissue detection from proteomics)
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class MatrisomeAdapter:
    def __init__(self, data_dir="template_package/data/matrisome"):
        self.data_dir = Path(data_dir)
        self.proteins = {}    # uniprot_id -> {gene, description}
        self.tissues = {}     # tissue_name -> set
        self.detections = []  # [{uniprot_id, tissue, species, reference}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load MatrisomeDB TSV files."""
        for fname in ['matrisome_human.tsv', 'matrisome_mouse.tsv']:
            path = self.data_dir / fname
            if path.exists():
                self._load_tsv(path)

        logger.info(f"MatrisomeDB: Loaded {len(self.proteins)} unique ECM proteins, "
                     f"{len(self.tissues)} tissues, {len(self.detections)} detections")

    def _load_tsv(self, path):
        """Load a MatrisomeDB TSV file."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene = row.get('gene', '').strip()
                    uniprot = row.get('uniprot', '').strip()
                    desc_raw = row.get('description', '').strip()
                    tissue = row.get('tissue', '').strip()
                    species = row.get('species', '').strip()
                    reference = row.get('reference', '').strip()

                    if not uniprot or not gene:
                        continue

                    # Clean description (remove HTML links)
                    desc = desc_raw.split('|')[0].strip() if '|' in desc_raw else desc_raw
                    desc = self._sanitize(desc)

                    # Register protein
                    if uniprot not in self.proteins:
                        self.proteins[uniprot] = {
                            'gene': gene,
                            'description': desc,
                            'species': species,
                        }

                    # Register tissue
                    if tissue:
                        self.tissues[tissue] = True

                    # Register detection
                    ref_clean = reference.split('|')[0].strip() if reference else ''
                    self.detections.append({
                        'uniprot_id': uniprot,
                        'tissue': tissue,
                        'species': species,
                        'reference': self._sanitize(ref_clean),
                    })
        except Exception as e:
            logger.warning(f"MatrisomeDB: Error loading {path}: {e}")

    def get_nodes(self):
        """
        Generate MatrisomeProtein nodes (unique ECM proteins).
        Yields: (id, label, properties)
        """
        logger.info("MatrisomeDB: Generating nodes...")
        count = 0

        for uniprot_id, data in self.proteins.items():
            props = {
                'gene_name': data['gene'],
                'description': data['description'],
                'species': data['species'],
                'source': 'MatrisomeDB',
            }

            yield (f"MATRISOME:{uniprot_id}", "MatrisomeProtein", props)
            count += 1

        logger.info(f"MatrisomeDB: Generated {count} MatrisomeProtein nodes")

    def get_edges(self):
        """
        Generate ProteinDetectedInTissue edges (deduplicated).
        Yields: (id, source, target, label, properties)
        """
        logger.info("MatrisomeDB: Generating edges...")
        seen = set()
        count = 0

        for det in self.detections:
            key = (det['uniprot_id'], det['tissue'])
            if key in seen:
                continue
            seen.add(key)

            if not det['tissue']:
                continue

            props = {
                'species': det['species'],
                'reference': det['reference'],
                'source': 'MatrisomeDB',
            }

            yield (
                None,
                f"MATRISOME:{det['uniprot_id']}",
                f"TISSUE:{det['tissue']}",
                "ProteinDetectedInTissue",
                props
            )
            count += 1

        logger.info(f"MatrisomeDB: Generated {count} ProteinDetectedInTissue edges")
