"""
miRTarBase Adapter for BioCypher.

Loads miRTarBase experimentally validated miRNA-target interactions and generates:
- MiRTarBaseInteraction edges (miRNA → target gene validated interactions)

miRTarBase is the largest curated database of experimentally validated
miRNA-target interactions, with evidence from reporter assays, Western blots, etc.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class MiRTarBaseAdapter:
    def __init__(self, data_dir="template_package/data/mirtarbase"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load miRTarBase validated interactions."""
        path = self.data_dir / 'hsa_MTI_v10.csv'
        if not path.exists():
            path = self.data_dir / 'hsa_MTI.csv'
        if not path.exists():
            logger.warning("miRTarBase: interaction data not found")
            return

        logger.info("miRTarBase: Loading miRNA-target interactions...")
        count = 0
        seen = set()

        with open(path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mirna = (row.get('miRNA') or '').strip()
                target_gene = (row.get('Target Gene') or '').strip()
                support_type = (row.get('Support Type') or '').strip()
                experiments = (row.get('Experiments') or '').strip()
                species_mirna = (row.get('Species (miRNA)') or '').strip()

                if not mirna or not target_gene:
                    continue
                if species_mirna != 'hsa':
                    continue

                # Deduplicate by miRNA-target pair, keeping strongest evidence
                key = (mirna, target_gene)
                if key in seen:
                    continue
                seen.add(key)

                self.interactions.append({
                    'mirna': mirna,
                    'target_gene': target_gene,
                    'support_type': support_type,
                    'experiments': experiments,
                })
                count += 1

        logger.info(f"miRTarBase: Loaded {count} unique miRNA-target interactions")

    def get_nodes(self):
        """No new nodes."""
        logger.info("miRTarBase: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate MiRTarBaseInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("miRTarBase: Generating edges...")
        count = 0

        for interaction in self.interactions:
            props = {
                'support_type': self._sanitize(interaction['support_type']),
                'experiments': self._sanitize(interaction['experiments']),
                'source': 'miRTarBase',
            }

            yield (
                None,
                interaction['mirna'],
                interaction['target_gene'],
                "MiRTarBaseInteraction",
                props
            )
            count += 1

        logger.info(f"miRTarBase: Generated {count} MiRTarBaseInteraction edges")
