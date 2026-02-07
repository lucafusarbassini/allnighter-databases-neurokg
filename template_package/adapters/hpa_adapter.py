"""
Human Protein Atlas (HPA) Adapter for BioCypher.

Loads HPA subcellular localization and tissue expression data and generates:
- SubcellularLocation nodes (cellular compartments)
- Tissue nodes (body tissues)
- GeneLocalizedTo edges (Gene → SubcellularLocation)
- GeneExpressedIn edges (Gene → Tissue)
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class HPAAdapter:
    def __init__(self, data_dir="template_package/data/hpa"):
        self.data_dir = data_dir
        self.gene_to_uniprot = {}  # gene_name -> uniprot_id
        self.locations = {}  # location_name -> GO_id
        self.tissues = set()
        self.localization_data = []
        self.expression_data = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load HPA data files."""
        self._load_gene_mapping()
        self._load_subcellular_location()
        self._load_tissue_expression()

    def _load_gene_mapping(self):
        """Load gene name → UniProt mapping from proteinatlas.tsv."""
        filepath = Path(self.data_dir) / 'proteinatlas.tsv'
        if not filepath.exists():
            logger.warning("HPA: proteinatlas.tsv not found")
            return

        logger.info("HPA: Loading gene-UniProt mappings...")
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('Gene', '').strip()
                uniprot = row.get('Uniprot', '').strip()
                if gene_name and uniprot:
                    self.gene_to_uniprot[gene_name] = uniprot

        logger.info(f"HPA: Mapped {len(self.gene_to_uniprot)} genes to UniProt IDs")

    def _load_subcellular_location(self):
        """Load subcellular location data."""
        filepath = Path(self.data_dir) / 'subcellular_location.tsv'
        if not filepath.exists():
            logger.warning("HPA: subcellular_location.tsv not found")
            return

        logger.info("HPA: Loading subcellular location data...")
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('Gene name', '').strip()
                reliability = row.get('Reliability', '').strip()
                main_loc = row.get('Main location', '').strip()
                additional_loc = row.get('Additional location', '').strip()
                go_id_str = row.get('GO id', '').strip()

                uniprot_id = self.gene_to_uniprot.get(gene_name)
                if not uniprot_id:
                    continue

                # Parse GO IDs from the GO id column
                # Format: "Cytosol (GO:0005829);Golgi apparatus (GO:0005794)"
                go_map = {}
                if go_id_str:
                    import re
                    for part in go_id_str.split(';'):
                        match = re.match(r'(.+?)\s*\((GO:\d+)\)', part.strip())
                        if match:
                            loc_name = match.group(1).strip()
                            go_id = match.group(2)
                            go_map[loc_name] = go_id
                            self.locations[loc_name] = go_id

                # Collect all locations
                all_locs = []
                if main_loc:
                    for loc in main_loc.split(';'):
                        loc = loc.strip()
                        if loc:
                            all_locs.append(('main', loc))
                if additional_loc:
                    for loc in additional_loc.split(';'):
                        loc = loc.strip()
                        if loc:
                            all_locs.append(('additional', loc))

                for loc_type, loc_name in all_locs:
                    self.localization_data.append({
                        'uniprot_id': uniprot_id,
                        'location': loc_name,
                        'go_id': go_map.get(loc_name, ''),
                        'location_type': loc_type,
                        'reliability': reliability,
                    })

        logger.info(f"HPA: Loaded {len(self.localization_data)} localization records, "
                     f"{len(self.locations)} unique locations")

    def _load_tissue_expression(self):
        """Load tissue expression data (protein level)."""
        filepath = Path(self.data_dir) / 'normal_tissue.tsv'
        if not filepath.exists():
            logger.warning("HPA: normal_tissue.tsv not found")
            return

        logger.info("HPA: Loading tissue expression data...")
        count = 0
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('Gene name', '').strip()
                tissue = row.get('Tissue', '').strip()
                cell_type = row.get('Cell type', '').strip()
                level = row.get('Level', '').strip()
                reliability = row.get('Reliability', '').strip()

                uniprot_id = self.gene_to_uniprot.get(gene_name)
                if not uniprot_id:
                    continue

                # Only keep detected proteins (Not detected = skip)
                if level == 'Not detected':
                    continue

                self.tissues.add(tissue)

                self.expression_data.append({
                    'uniprot_id': uniprot_id,
                    'tissue': tissue,
                    'cell_type': cell_type,
                    'level': level,
                    'reliability': reliability,
                })
                count += 1

        logger.info(f"HPA: Loaded {count} tissue expression records, "
                     f"{len(self.tissues)} unique tissues")

    def get_nodes(self):
        """
        Generate SubcellularLocation and Tissue nodes.
        Yields: (id, label, properties)
        """
        logger.info("HPA: Generating nodes...")
        loc_count = 0
        tissue_count = 0

        # SubcellularLocation nodes
        for loc_name, go_id in self.locations.items():
            node_id = go_id if go_id else f"HPA:{loc_name.replace(' ', '_')}"
            props = {
                'name': self._sanitize(loc_name),
                'go_id': go_id,
                'source': 'HPA',
            }
            yield (node_id, "SubcellularLocation", props)
            loc_count += 1

        # Tissue nodes
        for tissue in sorted(self.tissues):
            tissue_id = f"HPA_TISSUE:{tissue.replace(' ', '_')}"
            props = {
                'name': self._sanitize(tissue),
                'source': 'HPA',
            }
            yield (tissue_id, "Tissue", props)
            tissue_count += 1

        logger.info(f"HPA: Generated {loc_count} SubcellularLocation, "
                     f"{tissue_count} Tissue nodes")

    def get_edges(self):
        """
        Generate localization and expression edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("HPA: Generating edges...")
        loc_count = 0
        expr_count = 0

        # GeneLocalizedTo edges
        for record in self.localization_data:
            go_id = record['go_id']
            loc_id = go_id if go_id else f"HPA:{record['location'].replace(' ', '_')}"

            props = {
                'location_type': record['location_type'],
                'reliability': record['reliability'],
            }

            yield (
                None,
                record['uniprot_id'],
                loc_id,
                "GeneLocalizedTo",
                props
            )
            loc_count += 1

        # GeneExpressedIn edges (deduplicate by gene+tissue)
        seen = set()
        for record in self.expression_data:
            tissue_id = f"HPA_TISSUE:{record['tissue'].replace(' ', '_')}"
            key = (record['uniprot_id'], tissue_id)
            if key in seen:
                continue
            seen.add(key)

            props = {
                'level': record['level'],
                'reliability': record['reliability'],
                'source': 'HPA',
            }

            yield (
                None,
                record['uniprot_id'],
                tissue_id,
                "GeneExpressedIn",
                props
            )
            expr_count += 1

        logger.info(f"HPA: Generated {loc_count} GeneLocalizedTo, "
                     f"{expr_count} GeneExpressedIn edges")
