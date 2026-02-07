"""
Translocatome Adapter for BioCypher.

Loads Translocatome protein translocation data and generates:
- ProteinTranslocation edges (protein → translocation with confidence)

Translocatome catalogs proteins that translocate between subcellular
compartments, with localization scores and translocation confidence.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class TranslocatomeAdapter:
    def __init__(self, data_dir="template_package/data/translocatome"):
        self.data_dir = Path(data_dir)
        self.proteins = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Translocatome data."""
        path = self.data_dir / 'all_proteins.json'
        if not path.exists():
            logger.warning("Translocatome: data file not found")
            return

        logger.info("Translocatome: Loading translocation data...")

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        count = 0
        for entry in data:
            decision = entry.get('decision', '').strip()
            # Only include translocating proteins
            if 'non-translocating' in decision or decision == 'N/A':
                continue

            uniprot = entry.get('uniprotac', '').strip()
            gene_name = entry.get('gene_name', '').strip()
            localizations = entry.get('localizations', [])
            evidence_score = entry.get('evidence_score', 0)

            if not uniprot:
                continue

            # Parse localizations
            loc_str = '|'.join(str(l) for l in localizations[:5])

            self.proteins.append({
                'uniprot': uniprot,
                'gene_name': gene_name,
                'decision': decision,
                'localizations': loc_str,
                'evidence_score': evidence_score or 0,
            })
            count += 1

        logger.info(f"Translocatome: Loaded {count} translocating proteins")

    def get_nodes(self):
        """No new nodes."""
        logger.info("Translocatome: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate ProteinTranslocation edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("Translocatome: Generating edges...")
        count = 0

        for prot in self.proteins:
            props = {
                'gene_name': prot['gene_name'],
                'confidence': prot['decision'],
                'localizations': self._sanitize(prot['localizations']),
                'evidence_score': prot['evidence_score'],
                'source': 'Translocatome',
            }

            yield (
                None,
                prot['uniprot'],
                "TRANSLOCATION",
                "ProteinTranslocation",
                props
            )
            count += 1

        logger.info(f"Translocatome: Generated {count} ProteinTranslocation edges")
