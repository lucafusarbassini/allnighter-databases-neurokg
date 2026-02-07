"""
GOA Adapter for BioCypher.

Loads Gene Ontology Annotations (GOA) for human and generates:
- GOAnnotation edges (protein → GO term annotations)

GOA provides Gene Ontology annotations for UniProt proteins,
covering molecular function, biological process, and cellular component.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class GOAAdapter:
    def __init__(self, data_dir="template_package/data/goa"):
        self.data_dir = Path(data_dir)
        self.annotations = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load GOA human annotations."""
        path = self.data_dir / 'goa_human.gaf.gz'
        if not path.exists():
            logger.warning("GOA: annotation file not found")
            return

        logger.info("GOA: Loading human GO annotations...")
        count = 0
        seen = set()

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('!'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 15:
                    continue

                db = parts[0]
                uniprot_id = parts[1]
                symbol = parts[2]
                qualifier = parts[3]
                go_id = parts[4]
                evidence_code = parts[6]
                aspect = parts[8]  # F=function, P=process, C=component

                if db != 'UniProtKB':
                    continue

                # Deduplicate by protein-GO pair
                key = (uniprot_id, go_id)
                if key in seen:
                    continue
                seen.add(key)

                self.annotations.append({
                    'uniprot_id': uniprot_id,
                    'symbol': symbol,
                    'qualifier': qualifier,
                    'go_id': go_id,
                    'evidence_code': evidence_code,
                    'aspect': aspect,
                })
                count += 1

        logger.info(f"GOA: Loaded {count} unique protein-GO annotations")

    def get_nodes(self):
        """No new nodes."""
        logger.info("GOA: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate GOAnnotation edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("GOA: Generating edges...")
        count = 0

        for ann in self.annotations:
            props = {
                'gene_symbol': self._sanitize(ann['symbol']),
                'qualifier': self._sanitize(ann['qualifier']),
                'evidence_code': ann['evidence_code'],
                'aspect': ann['aspect'],
                'source': 'GOA',
            }

            yield (
                None,
                ann['uniprot_id'],
                ann['go_id'],
                "GOAnnotation",
                props
            )
            count += 1

        logger.info(f"GOA: Generated {count} GOAnnotation edges")
