"""
REPAIRtoire Adapter for BioCypher.

Loads REPAIRtoire DNA repair protein data and generates:
- DNARepairProtein nodes (repair proteins with pathway assignments)
- ProteinInRepairPathway edges (protein → repair pathway)

REPAIRtoire catalogs proteins involved in DNA damage repair pathways
including BER, NER, MMR, HR, NHEJ, DRR, TLS, and FA.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class RepairToireAdapter:
    def __init__(self, data_dir="template_package/data/repairtoire"):
        self.data_dir = Path(data_dir)
        self.proteins = []
        self.pathways = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load REPAIRtoire human repair proteins."""
        path = self.data_dir / 'human_repair_proteins.json'
        if not path.exists():
            logger.warning("REPAIRtoire: data file not found")
            return

        logger.info("REPAIRtoire: Loading DNA repair protein data...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        count = 0
        for entry in data:
            name = entry.get('name', '').strip()
            full_name = entry.get('full_name', '').strip()
            organism = entry.get('organism', '').strip()
            uniprot_ids = entry.get('uniprot_ids', [])
            pathways = entry.get('pathways', [])
            function_desc = entry.get('function', '').strip()

            if not name or not uniprot_ids:
                continue

            # Collect pathway names
            pathway_names = []
            for pw in pathways:
                pw_name = pw.get('name', '').strip()
                pw_id = pw.get('id', '').strip()
                if pw_name:
                    pathway_names.append(pw_name)
                    self.pathways[pw_id] = pw_name

            self.proteins.append({
                'name': name,
                'full_name': full_name,
                'uniprot_ids': uniprot_ids,
                'pathways': pathway_names,
                'function': function_desc,
            })
            count += 1

        logger.info(f"REPAIRtoire: Loaded {count} DNA repair proteins in {len(self.pathways)} pathways")

    def get_nodes(self):
        """No new nodes - uses existing Gene nodes via UniProt ID."""
        logger.info("REPAIRtoire: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate ProteinInRepairPathway edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("REPAIRtoire: Generating edges...")
        count = 0

        for protein in self.proteins:
            for uniprot_id in protein['uniprot_ids']:
                for pathway_name in protein['pathways']:
                    props = {
                        'gene_name': protein['name'],
                        'pathway_name': pathway_name,
                        'function': self._sanitize(protein['function'][:500]),
                        'source': 'REPAIRtoire',
                    }

                    yield (
                        None,
                        uniprot_id,
                        f"REPAIR_PATHWAY:{pathway_name}",
                        "ProteinInRepairPathway",
                        props
                    )
                    count += 1

        logger.info(f"REPAIRtoire: Generated {count} ProteinInRepairPathway edges")
