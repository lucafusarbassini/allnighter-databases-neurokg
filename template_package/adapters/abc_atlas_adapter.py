"""
Allen Brain Cell (ABC) Atlas Adapter for BioCypher.

Parses Allen Brain Atlas API JSON files and cell type CSV data.
Generates:
- BrainRegion nodes (from allen_brain_structures.json)
- CellType nodes (from allen_mouse_cell_types.json, p2, and CSV)
- Gene nodes (from allen_brain_genes.json)

Data files:
- allen_brain_structures.json (1,327 brain regions)
- allen_mouse_cell_types.json + p2 (2,333 cell type specimens)
- allen_brain_genes.json (2,000 genes)
- allen_cell_types_database.csv (73,348 cell entries)
"""

import csv
import json
from pathlib import Path
from biocypher._logger import logger


class ABCAtlasAdapter:
    def __init__(self, data_dir="template_package/data/abc_atlas"):
        self.data_dir = Path(data_dir)
        self.brain_structures = []
        self.cell_type_specimens = []
        self.genes = []
        self.csv_cell_types = {}  # deduplicated by cell_type_accession_id
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_json_api(self, filename):
        """Load Allen API JSON (has 'msg' key with list of records)."""
        path = self.data_dir / filename
        if not path.exists():
            logger.warning(f"ABC Atlas: {filename} not found")
            return []
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if isinstance(data, dict) and 'msg' in data:
            return data['msg']
        elif isinstance(data, list):
            return data
        return []

    def _load_data(self):
        """Load all Allen Brain Atlas data files."""
        # Brain structures
        self.brain_structures = self._load_json_api('allen_brain_structures.json')
        logger.info(f"ABC Atlas: Loaded {len(self.brain_structures)} brain structures")

        # Cell type specimens (two pages)
        p1 = self._load_json_api('allen_mouse_cell_types.json')
        p2 = self._load_json_api('allen_mouse_cell_types_p2.json')
        self.cell_type_specimens = p1 + p2
        logger.info(f"ABC Atlas: Loaded {len(self.cell_type_specimens)} cell type specimens")

        # Genes
        self.genes = self._load_json_api('allen_brain_genes.json')
        logger.info(f"ABC Atlas: Loaded {len(self.genes)} genes")

        # CSV cell types (deduplicate by cell_type_accession_id)
        csv_path = self.data_dir / 'allen_cell_types_database.csv'
        if csv_path.exists():
            with open(csv_path, 'r', encoding='utf-8', errors='replace') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    ct_id = row.get('cell_type_accession_id', '').strip()
                    if ct_id and ct_id not in self.csv_cell_types:
                        self.csv_cell_types[ct_id] = {
                            'cell_type_accession_id': ct_id,
                            'cluster_label': row.get('cluster_label', '').strip(),
                            'class_label': row.get('class_label', '').strip(),
                            'subclass_label': row.get('subclass_label', '').strip(),
                            'neighborhood_label': row.get('neighborhood_label', '').strip(),
                            'region_label': row.get('region_label', '').strip(),
                            'cell_type_alias_label': row.get('cell_type_alias_label', '').strip(),
                            'cell_type_designation_label': row.get('cell_type_designation_label', '').strip(),
                            'cluster_color': row.get('cluster_color', '').strip(),
                        }
            logger.info(f"ABC Atlas: Loaded {len(self.csv_cell_types)} unique cell types from CSV")

    def get_nodes(self):
        """
        Yield BrainRegion nodes, CellType nodes, and Gene nodes.
        Each: (id, label, properties)
        """
        # --- BrainRegion nodes ---
        logger.info("ABC Atlas: Generating BrainRegion nodes...")
        br_count = 0
        for s in self.brain_structures:
            sid = s.get('id')
            if sid is None:
                continue
            props = {
                'name': self._sanitize(s.get('name', '')),
                'acronym': self._sanitize(s.get('acronym', '')),
                'color_hex': self._sanitize(s.get('color_hex_triplet', '')),
                'graph_order': s.get('graph_order', 0),
                'structure_level': s.get('st_level', 0),
                'depth': s.get('depth', 0),
                'parent_structure_id': s.get('parent_structure_id'),
                'structure_id_path': self._sanitize(s.get('structure_id_path', '')),
                'ontology_id': s.get('ontology_id', 0),
                'hemisphere_id': s.get('hemisphere_id', 3),
                'source': 'Allen Brain Atlas',
            }
            yield (f"allen:structure:{sid}", "BrainRegion", props)
            br_count += 1
        logger.info(f"ABC Atlas: Generated {br_count} BrainRegion nodes")

        # --- CellType nodes from CSV (deduplicated) ---
        logger.info("ABC Atlas: Generating CellType nodes from CSV...")
        ct_count = 0
        for ct_id, ct in self.csv_cell_types.items():
            props = {
                'name': self._sanitize(ct['cluster_label']),
                'class_label': self._sanitize(ct['class_label']),
                'subclass_label': self._sanitize(ct['subclass_label']),
                'neighborhood_label': self._sanitize(ct['neighborhood_label']),
                'region': self._sanitize(ct['region_label']),
                'alias': self._sanitize(ct['cell_type_alias_label']),
                'designation': self._sanitize(ct['cell_type_designation_label']),
                'color_hex': self._sanitize(ct['cluster_color']),
                'source': 'Allen Cell Types Database',
            }
            yield (f"abc:{ct_id}", "CellType", props)
            ct_count += 1
        logger.info(f"ABC Atlas: Generated {ct_count} CellType nodes from CSV")

        # --- Gene nodes ---
        logger.info("ABC Atlas: Generating Gene nodes...")
        gene_count = 0
        for g in self.genes:
            gid = g.get('id')
            if gid is None:
                continue
            props = {
                'symbol': self._sanitize(g.get('acronym', '')),
                'name': self._sanitize(g.get('name', '')),
                'entrez_id': g.get('entrez_id'),
                'ensembl_id': self._sanitize(g.get('ensembl_id', '') or ''),
                'chromosome_id': g.get('chromosome_id'),
                'organism_id': g.get('organism_id'),
                'version_status': self._sanitize(g.get('version_status', '')),
                'source': 'Allen Brain Atlas',
            }
            yield (f"allen:gene:{gid}", "Gene", props)
            gene_count += 1
        logger.info(f"ABC Atlas: Generated {gene_count} Gene nodes")

    def get_edges(self):
        """
        Yield CellTypeSpecimen edges linking specimens to donor/species context.
        Each: (id, source, target, label, properties)
        """
        logger.info("ABC Atlas: Generating CellTypeSpecimen edges...")
        count = 0
        for spec in self.cell_type_specimens:
            spec_id = spec.get('specimen__id')
            donor_id = spec.get('donor__id')
            if not spec_id or not donor_id:
                continue
            species = self._sanitize(spec.get('donor__species', ''))
            props = {
                'donor_name': self._sanitize(spec.get('donor__name', '')),
                'donor_sex': self._sanitize(spec.get('donor__sex', '')),
                'donor_age': self._sanitize(spec.get('donor__age', '')),
                'species': species,
                'disease_state': self._sanitize(spec.get('donor__disease_state', '')),
                'cell_reporter_status': self._sanitize(spec.get('cell_reporter_status', '') or ''),
                'source': 'Allen Cell Types',
            }
            yield (
                f"allen:specimen:{spec_id}",
                f"allen:specimen:{spec_id}",
                f"allen:donor:{donor_id}",
                "CellTypeSpecimen",
                props,
            )
            count += 1
        logger.info(f"ABC Atlas: Generated {count} CellTypeSpecimen edges")
