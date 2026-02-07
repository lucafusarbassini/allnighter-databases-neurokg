"""
GPCRdb (G Protein-Coupled Receptor Database) Adapter for BioCypher.

Loads GPCR data from GPCRdb JSON API dumps and generates:
- GPCRFamily nodes (receptor classes and families)
- GPCRLigandInteraction edges (Gene → ligand name)
"""

import json
import re
from pathlib import Path
from biocypher._logger import logger


class GPCRdbAdapter:
    def __init__(self, data_dir="template_package/data/gpcrdb"):
        self.data_dir = data_dir
        self.receptors = []
        self.families = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        # Remove HTML tags
        text = re.sub(r'<[^>]+>', '', text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load GPCRdb JSON data."""
        # Load human receptors
        receptors_path = Path(self.data_dir) / 'human_receptors.json'
        if receptors_path.exists():
            logger.info("GPCRdb: Loading human receptors...")
            with open(receptors_path, 'r') as f:
                self.receptors = json.load(f)
            logger.info(f"GPCRdb: Loaded {len(self.receptors)} human GPCRs")

        # Load protein families
        families_path = Path(self.data_dir) / 'protein_families.json'
        if families_path.exists():
            logger.info("GPCRdb: Loading protein families...")
            with open(families_path, 'r') as f:
                families_data = json.load(f)
                # Build family hierarchy
                self._parse_families(families_data)
            logger.info(f"GPCRdb: Loaded {len(self.families)} families")

    def _parse_families(self, data, prefix=""):
        """Recursively parse protein family hierarchy."""
        if isinstance(data, list):
            for item in data:
                self._parse_families(item, prefix)
        elif isinstance(data, dict):
            for key, value in data.items():
                family_id = f"GPCR:{key}" if not key.startswith("GPCR:") else key
                self.families[family_id] = {
                    'name': self._sanitize(key),
                    'parent': prefix if prefix else None,
                }
                if isinstance(value, dict):
                    self._parse_families(value, family_id)
                elif isinstance(value, list):
                    self._parse_families(value, family_id)

    def get_nodes(self):
        """
        Generate GPCR-related nodes.
        Yields: (id, label, properties)
        """
        logger.info("GPCRdb: Generating nodes...")
        count = 0

        # GPCR Family nodes
        seen_classes = set()
        for receptor in self.receptors:
            receptor_class = receptor.get('receptor_class', '')
            receptor_family = receptor.get('receptor_family', '')
            ligand_type = receptor.get('ligand_type', '')

            for category, cat_type in [
                (receptor_class, 'class'),
                (receptor_family, 'family'),
                (ligand_type, 'ligand_type'),
            ]:
                if category and category not in seen_classes:
                    seen_classes.add(category)
                    cat_id = f"GPCR:{self._sanitize(category).replace(' ', '_')}"
                    props = {
                        'name': self._sanitize(category),
                        'category_type': cat_type,
                        'source': 'GPCRdb',
                    }
                    yield (cat_id, "GPCRFamily", props)
                    count += 1

        logger.info(f"GPCRdb: Generated {count} GPCRFamily nodes")

    def get_edges(self):
        """
        Generate GPCR-related edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("GPCRdb: Generating edges...")
        classify_count = 0
        ligand_count = 0

        for receptor in self.receptors:
            accession = receptor.get('accession', '')
            if not accession:
                continue

            receptor_class = receptor.get('receptor_class', '')
            receptor_family = receptor.get('receptor_family', '')

            # Classification edges (Gene → GPCRFamily)
            if receptor_class:
                class_id = f"GPCR:{self._sanitize(receptor_class).replace(' ', '_')}"
                yield (
                    None,
                    accession,
                    class_id,
                    "GPCRClassifiedAs",
                    {'level': 'class'}
                )
                classify_count += 1

            if receptor_family:
                family_id = f"GPCR:{self._sanitize(receptor_family).replace(' ', '_')}"
                yield (
                    None,
                    accession,
                    family_id,
                    "GPCRClassifiedAs",
                    {'level': 'family'}
                )
                classify_count += 1

            # Endogenous ligand edges
            ligands = receptor.get('endogenous_ligands', [])
            if isinstance(ligands, list):
                for ligand in ligands:
                    if isinstance(ligand, dict):
                        ligand_name = self._sanitize(ligand.get('name', ''))
                        if ligand_name:
                            yield (
                                None,
                                accession,
                                f"LIGAND:{ligand_name.replace(' ', '_')}",
                                "GPCRBindsLigand",
                                {'ligand_name': ligand_name}
                            )
                            ligand_count += 1

        logger.info(f"GPCRdb: Generated {classify_count} GPCRClassifiedAs, "
                     f"{ligand_count} GPCRBindsLigand edges")
