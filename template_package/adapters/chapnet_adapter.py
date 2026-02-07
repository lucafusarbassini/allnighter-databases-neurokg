"""
ChapNet Adapter for BioCypher.

Loads chaperone-protein co-expression/interaction network data from ChapNet
(https://netbio.bgu.ac.il/chapnet/). Data is in Cytoscape JSON format
containing chaperone correlation networks across tissues.

Generates:
- Chaperone interaction edges (protein-protein co-expression correlations)
"""

import json
import csv
from pathlib import Path
from biocypher._logger import logger


class ChapNetAdapter:
    def __init__(self, data_dir="template_package/data/chapnet"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self.genes = set()
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ChapNet correlation data."""
        # Try JSON files first (Cytoscape format)
        json_files = sorted(self.data_dir.glob('*.json'))
        if json_files:
            self._load_json_networks(json_files)

        # Also load CSV correlation data
        corr_path = self.data_dir / 'ChaperoneCorrelation.csv'
        if corr_path.exists():
            self._load_correlation_csv(corr_path)

    def _load_json_networks(self, json_files):
        """Load Cytoscape JSON network files."""
        for jf in json_files:
            try:
                with open(jf, 'r', encoding='utf-8') as f:
                    data = json.load(f)

                elements = data.get('elements', {})
                edges = elements.get('edges', [])
                nodes = elements.get('nodes', [])

                # Track gene Ensembl IDs from nodes
                for node in nodes:
                    nd = node.get('data', {})
                    ensembl = nd.get('Ensembl', '')
                    if ensembl:
                        self.genes.add(ensembl)

                network_name = jf.stem

                for edge in edges:
                    ed = edge.get('data', {})
                    source = ed.get('source', '')
                    target = ed.get('target', '')
                    source_sym = ed.get('sourceSymbol', ed.get('Protein1Symbol', '')).strip()
                    target_sym = ed.get('targetSymbol', ed.get('Protein2Symbol', '')).strip()

                    if not source or not target:
                        continue

                    self.interactions.append({
                        'source': source,
                        'target': target,
                        'source_symbol': source_sym,
                        'target_symbol': target_sym,
                        'network': network_name,
                    })

                logger.info(f"ChapNet: Loaded {len(edges)} edges from {jf.name}")
            except Exception as e:
                logger.warning(f"ChapNet: Error loading {jf.name}: {e}")

    def _load_correlation_csv(self, path):
        """Load bulk correlation CSV (tissue-specific correlations)."""
        try:
            count = 0
            seen = set()
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    p1 = row.get('Protein1', '').strip()
                    p2 = row.get('Protein2', '').strip()
                    p1_sym = row.get('Protein1Symbol', '').strip()
                    p2_sym = row.get('Protein2Symbol', '').strip()

                    if not p1 or not p2:
                        continue

                    pair = tuple(sorted([p1, p2]))
                    if pair in seen:
                        continue
                    seen.add(pair)

                    self.interactions.append({
                        'source': p1,
                        'target': p2,
                        'source_symbol': p1_sym,
                        'target_symbol': p2_sym,
                        'network': 'ChaperoneCorrelation',
                    })
                    count += 1

            logger.info(f"ChapNet: Loaded {count} unique pairs from correlation CSV")
        except Exception as e:
            logger.warning(f"ChapNet: Error loading correlation CSV: {e}")

    def get_nodes(self):
        """No additional nodes - uses existing Gene nodes."""
        logger.info("ChapNet: No additional nodes to generate (uses Gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate chaperone interaction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info(f"ChapNet: Generating edges from {len(self.interactions)} interactions...")
        count = 0
        seen = set()

        for inter in self.interactions:
            source = inter['source']
            target = inter['target']
            pair = tuple(sorted([source, target]))
            if pair in seen:
                continue
            seen.add(pair)

            props = {
                'source_symbol': self._sanitize(inter['source_symbol']),
                'target_symbol': self._sanitize(inter['target_symbol']),
                'network': inter['network'],
                'source_db': 'ChapNet',
            }

            yield (None, source, target, "ChaperoneInteraction", props)
            count += 1

        logger.info(f"ChapNet: Generated {count} chaperone interaction edges")
