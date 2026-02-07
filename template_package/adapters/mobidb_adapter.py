"""
MobiDB Adapter for BioCypher.

Loads MobiDB intrinsic disorder annotations and generates:
- DisorderRegion edges (protein → disorder annotation with region coordinates)

MobiDB integrates curated and predicted intrinsically disordered regions (IDRs)
from DisProt, IDEAL, and prediction tools.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class MobiDBAdapter:
    def __init__(self, data_dir="template_package/data/mobidb"):
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
        """Load MobiDB compact JSONL data."""
        path = self.data_dir / 'mobidb_human_disorder_compact.jsonl'
        if not path.exists():
            logger.warning("MobiDB: data file not found")
            return

        logger.info("MobiDB: Loading intrinsic disorder data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:
                    continue

                acc = entry.get('acc', '')
                gene = entry.get('gene', '')
                name = entry.get('name', '')
                length = entry.get('length', 0)

                # Extract curated disorder regions (priority consensus)
                disorder_info = entry.get('curated-disorder-priority', {})
                if not disorder_info:
                    # Fall back to predicted
                    disorder_info = entry.get('prediction-disorder-th_50', {})

                regions = disorder_info.get('regions', [])
                content_fraction = disorder_info.get('content_fraction', 0.0)
                content_count = disorder_info.get('content_count', 0)

                if not acc or not regions:
                    continue

                self.proteins.append({
                    'acc': acc,
                    'gene': gene,
                    'name': name,
                    'length': length,
                    'regions': regions,
                    'content_fraction': content_fraction,
                    'content_count': content_count,
                })
                count += 1

        logger.info(f"MobiDB: Loaded {count} proteins with disorder regions")

    def get_nodes(self):
        """No new nodes - links to existing Gene nodes via UniProt ID."""
        logger.info("MobiDB: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate ProteinHasDisorder edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("MobiDB: Generating edges...")
        count = 0

        for protein in self.proteins:
            # Create one edge per protein summarizing its disorder
            regions_str = '|'.join(f"{r[0]}-{r[1]}" for r in protein['regions'])

            props = {
                'disorder_regions': regions_str,
                'num_regions': len(protein['regions']),
                'content_fraction': protein['content_fraction'],
                'disordered_residues': protein['content_count'],
                'protein_length': protein['length'],
                'source': 'MobiDB',
            }

            yield (
                None,
                protein['acc'],
                f"IDR:{protein['acc']}",
                "ProteinHasDisorder",
                props
            )
            count += 1

        logger.info(f"MobiDB: Generated {count} ProteinHasDisorder edges")
