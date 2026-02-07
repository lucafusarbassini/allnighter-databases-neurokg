"""
ELM (Eukaryotic Linear Motif) Adapter for BioCypher.

Loads ELM motif classes and instances, generating:
- LinearMotif nodes (motif classes with regex patterns)
- ProteinHasMotif edges (protein → motif instance mapping)
"""

import csv
import re
from pathlib import Path
from biocypher._logger import logger


class ELMAdapter:
    def __init__(self, data_dir="template_package/data/elm"):
        self.data_dir = Path(data_dir)
        self.classes = []
        self.instances = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = re.sub(r'<[^>]+>', '', text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ELM TSV files."""
        # Load ELM classes (motif definitions)
        classes_path = self.data_dir / 'elms.tsv'
        if classes_path.exists():
            self._load_tsv(classes_path, self.classes, 'classes')

        # Fallback: try elm_classes.tsv if elms.tsv was HTML
        if not self.classes:
            classes_path2 = self.data_dir / 'elm_classes.tsv'
            if classes_path2.exists():
                self._load_tsv(classes_path2, self.classes, 'classes')

        # Load ELM instances
        instances_path = self.data_dir / 'elm_instances.tsv'
        if instances_path.exists():
            self._load_tsv(instances_path, self.instances, 'instances')

        logger.info(f"ELM: Loaded {len(self.classes)} motif classes, {len(self.instances)} instances")

    def _load_tsv(self, path, target_list, data_type):
        """Load a TSV file, skipping comment lines."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                # Skip comment lines
                lines = []
                header = None
                for line in f:
                    line = line.strip()
                    if line.startswith('#'):
                        continue
                    if not line:
                        continue
                    # Check if this is HTML (failed download)
                    if line.startswith('<!DOCTYPE') or line.startswith('<html'):
                        logger.warning(f"ELM: {path.name} appears to be HTML, skipping")
                        return
                    if header is None:
                        header = line
                        lines.append(header)
                    else:
                        lines.append(line)

                if not lines:
                    return

                reader = csv.DictReader(lines, delimiter='\t', quotechar='"')
                for row in reader:
                    target_list.append(row)

                logger.info(f"ELM: Loaded {len(target_list)} {data_type} from {path.name}")
        except Exception as e:
            logger.warning(f"ELM: Error loading {path}: {e}")

    def get_nodes(self):
        """
        Generate ELM linear motif nodes.
        Yields: (id, label, properties)
        """
        logger.info("ELM: Generating nodes...")
        count = 0

        for cls in self.classes:
            acc = cls.get('Accession', cls.get('"Accession"', ''))
            elm_id = cls.get('ELMIdentifier', cls.get('"ELMIdentifier"', ''))
            func_name = cls.get('FunctionalSiteName', cls.get('"FunctionalSiteName"', ''))
            description = cls.get('Description', cls.get('"Description"', ''))
            regex = cls.get('Regex', cls.get('"Regex"', ''))
            probability = cls.get('Probability', cls.get('"Probability"', ''))
            num_instances = cls.get('#Instances', cls.get('"#Instances"', ''))
            num_pdb = cls.get('#Instances_in_PDB', cls.get('"#Instances_in_PDB"', ''))

            if not acc:
                continue

            node_id = f"ELM:{acc}"

            # Determine motif type from identifier prefix
            motif_type = ''
            if elm_id:
                parts = elm_id.split('_')
                if parts:
                    motif_type = parts[0]  # CLV, DEG, DOC, LIG, MOD, TRG

            props = {
                'elm_identifier': self._sanitize(elm_id),
                'name': self._sanitize(func_name),
                'description': self._sanitize(description),
                'regex': self._sanitize(regex),
                'motif_type': motif_type,
                'probability': float(probability) if probability else 0.0,
                'num_instances': int(num_instances) if num_instances else 0,
                'num_pdb_instances': int(num_pdb) if num_pdb else 0,
                'source': 'ELM',
            }

            yield (node_id, "LinearMotif", props)
            count += 1

        logger.info(f"ELM: Generated {count} LinearMotif nodes")

    def get_edges(self):
        """
        Generate ELM edges linking proteins to motif instances.
        Yields: (id, source, target, label, properties)
        """
        logger.info("ELM: Generating edges...")
        count = 0

        for inst in self.instances:
            primary_acc = inst.get('Primary_Acc', inst.get('"Primary_Acc"', ''))
            elm_identifier = inst.get('ELMIdentifier', inst.get('"ELMIdentifier"', ''))
            elm_acc = inst.get('Accession', inst.get('"Accession"', ''))
            organism = inst.get('Organism', inst.get('"Organism"', ''))
            start = inst.get('Start', inst.get('"Start"', ''))
            end = inst.get('End', inst.get('"End"', ''))
            logic = inst.get('InstanceLogic', inst.get('"InstanceLogic"', ''))
            methods = inst.get('Methods', inst.get('"Methods"', ''))

            if not primary_acc or not elm_identifier:
                continue

            # Filter to human and mouse
            org_clean = self._sanitize(organism)
            if 'Homo sapiens' not in org_clean and 'Mus musculus' not in org_clean:
                continue

            # Find the class accession for this identifier
            class_acc = None
            for cls in self.classes:
                cls_id = cls.get('ELMIdentifier', cls.get('"ELMIdentifier"', ''))
                cls_acc = cls.get('Accession', cls.get('"Accession"', ''))
                if cls_id == elm_identifier:
                    class_acc = cls_acc
                    break

            target_id = f"ELM:{class_acc}" if class_acc else f"ELM:{elm_identifier}"

            props = {
                'instance_accession': self._sanitize(elm_acc),
                'start_position': int(start) if start else 0,
                'end_position': int(end) if end else 0,
                'instance_logic': self._sanitize(logic),
                'detection_methods': self._sanitize(methods),
                'organism': org_clean,
            }

            yield (None, primary_acc, target_id, "ProteinHasMotif", props)
            count += 1

        logger.info(f"ELM: Generated {count} ProteinHasMotif edges")
