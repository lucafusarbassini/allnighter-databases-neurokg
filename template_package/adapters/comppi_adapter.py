"""
ComPPI (Compartmentalized Protein-Protein Interaction) Adapter for BioCypher.

Loads ComPPI human protein localization data and generates:
- CompartmentalizedInteraction edges (scored PPI with compartment context)

ComPPI integrates PPIs with subcellular localization data from multiple
databases, assigning interaction scores based on compartmental evidence.
"""

import gzip
import csv
from pathlib import Path
from biocypher._logger import logger


class ComPPIAdapter:
    def __init__(self, data_dir="template_package/data/comppi"):
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
        """Load ComPPI integrated PPI data."""
        ppi_path = self.data_dir / 'comppi_integrated_ppi_hsapiens.txt.gz'
        if not ppi_path.exists():
            logger.warning("ComPPI: PPI file not found")
            return

        logger.info("ComPPI: Loading compartmentalized PPIs...")
        count = 0
        skipped = 0

        with gzip.open(ppi_path, 'rt', encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                prot_a = row.get('Protein A', '').strip()
                prot_b = row.get('Protein B', '').strip()
                score = row.get('Interaction Score', '').strip()
                sys_type = row.get('Interaction Experimental System Type', '').strip()
                source_db = row.get('Interaction Source Database', '').strip()

                if not prot_a or not prot_b:
                    skipped += 1
                    continue

                try:
                    score_val = float(score)
                except (ValueError, TypeError):
                    score_val = 0.0

                self.interactions.append({
                    'prot_a': prot_a,
                    'prot_b': prot_b,
                    'score': score_val,
                    'system_type': sys_type,
                    'source_db': source_db,
                })
                count += 1

                # Cap at 500K for memory
                if count >= 500000:
                    break

        logger.info(f"ComPPI: Loaded {count} interactions (skipped {skipped})")

    def get_nodes(self):
        """No new nodes - uses existing Gene nodes."""
        logger.info("ComPPI: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate CompartmentalizedInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("ComPPI: Generating edges...")
        seen = set()
        count = 0

        for inter in self.interactions:
            # Deduplicate A-B / B-A
            key = tuple(sorted([inter['prot_a'], inter['prot_b']]))
            if key in seen:
                continue
            seen.add(key)

            props = {
                'interaction_score': inter['score'],
                'system_type': self._sanitize(inter['system_type']),
                'source_db': self._sanitize(inter['source_db']),
                'source': 'ComPPI',
            }

            yield (
                None,
                inter['prot_a'],
                inter['prot_b'],
                "CompartmentalizedInteraction",
                props
            )
            count += 1

        logger.info(f"ComPPI: Generated {count} CompartmentalizedInteraction edges")
