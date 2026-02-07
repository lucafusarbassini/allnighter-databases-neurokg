"""
DGIdb (Drug Gene Interaction Database) Adapter for BioCypher.

Loads DGIdb drug-gene interaction data from GraphQL JSON exports and generates:
- Drug nodes (approved and investigational drugs)
- DrugGeneInteraction edges (drug → gene with interaction type and score)
"""

import json
from pathlib import Path
from biocypher._logger import logger


class DGIdbAdapter:
    def __init__(self, data_dir="template_package/data/dgidb"):
        self.data_dir = Path(data_dir)
        self.drugs = {}           # concept_id -> {name, approved}
        self.interactions = []    # [{drug_id, gene_name, score, types, sources}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load DGIdb JSON files."""
        interactions_path = self.data_dir / 'interactions.json'
        if not interactions_path.exists():
            logger.warning("DGIdb: interactions.json not found")
            return

        logger.info("DGIdb: Loading drug-gene interactions...")

        with open(interactions_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        genes = data.get('data', {}).get('genes', {}).get('nodes', [])

        for gene in genes:
            gene_name = gene.get('name', '')
            if not gene_name:
                continue

            for inter in gene.get('interactions', []):
                drug = inter.get('drug', {})
                drug_name = drug.get('name', '')
                drug_id = drug.get('conceptId', '')
                approved = drug.get('approved', False)

                if not drug_name:
                    continue

                # Use concept ID as primary ID, fall back to name
                if drug_id:
                    node_id = drug_id
                else:
                    node_id = f"DGIDB:{drug_name}"

                # Register drug
                if node_id not in self.drugs:
                    self.drugs[node_id] = {
                        'name': drug_name,
                        'approved': approved,
                    }

                # Extract interaction metadata
                score = inter.get('interactionScore', 0.0)
                types = inter.get('interactionTypes', [])
                type_str = '|'.join(t.get('type', '') for t in types if t.get('type'))
                directionality = ''
                if types:
                    directionality = types[0].get('directionality', '')

                sources = inter.get('sources', [])
                source_str = '|'.join(s.get('fullName', '') for s in sources if s.get('fullName'))

                self.interactions.append({
                    'drug_id': node_id,
                    'gene_name': gene_name,
                    'score': score,
                    'interaction_types': type_str,
                    'directionality': directionality,
                    'sources': source_str,
                })

        logger.info(f"DGIdb: Loaded {len(self.drugs)} drugs, {len(self.interactions)} interactions")

    def get_nodes(self):
        """
        Generate Drug nodes.
        Yields: (id, label, properties)
        """
        logger.info("DGIdb: Generating drug nodes...")
        count = 0

        for drug_id, data in self.drugs.items():
            props = {
                'name': self._sanitize(data['name']),
                'approved': data.get('approved', False),
                'source': 'DGIdb',
            }

            yield (drug_id, "Drug", props)
            count += 1

        logger.info(f"DGIdb: Generated {count} Drug nodes")

    def get_edges(self):
        """
        Generate DrugGeneInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("DGIdb: Generating edges...")
        count = 0

        for inter in self.interactions:
            props = {
                'interaction_score': inter['score'],
                'interaction_types': self._sanitize(inter['interaction_types']),
                'directionality': inter['directionality'],
                'sources': self._sanitize(inter['sources']),
            }

            yield (
                None,
                inter['drug_id'],
                inter['gene_name'],
                "DrugGeneInteraction",
                props
            )
            count += 1

        logger.info(f"DGIdb: Generated {count} DrugGeneInteraction edges")
