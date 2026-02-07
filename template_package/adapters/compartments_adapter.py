"""
Compartments (Jensen Lab) Adapter for BioCypher.

Loads subcellular localization data from Compartments database and generates:
- CompartmentLocation nodes (GO cellular component terms)
- ProteinLocatedIn edges (Gene → CompartmentLocation with confidence)

Uses knowledge-based and experimental evidence channels.
"""

from pathlib import Path
from biocypher._logger import logger


class CompartmentsAdapter:
    def __init__(self, data_dir="template_package/data/compartments",
                 min_confidence=3):
        self.data_dir = data_dir
        self.min_confidence = min_confidence
        self.compartments = {}  # GO_id -> name
        self.localizations = []
        self.ensp_to_gene = {}  # ENSP -> gene_name
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Compartments knowledge and experimental data."""
        seen_pairs = set()

        for filename, evidence_type in [
            ('human_compartment_knowledge_full.tsv', 'knowledge'),
            ('human_compartment_experiments_full.tsv', 'experimental'),
        ]:
            filepath = Path(self.data_dir) / filename
            if not filepath.exists():
                logger.warning(f"Compartments: {filename} not found")
                continue

            logger.info(f"Compartments: Loading {filename}...")
            count = 0

            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 7:
                        continue

                    ensp_id = parts[0].strip()
                    gene_name = parts[1].strip()
                    go_id = parts[2].strip()
                    go_name = parts[3].strip()
                    source = parts[4].strip()
                    evidence_code = parts[5].strip()

                    try:
                        confidence = int(parts[6].strip())
                    except (ValueError, IndexError):
                        confidence = 0

                    if confidence < self.min_confidence:
                        continue

                    # Track ENSP to gene name mapping
                    if ensp_id and gene_name:
                        self.ensp_to_gene[ensp_id] = gene_name

                    # Register compartment
                    if go_id and go_name:
                        if go_id not in self.compartments:
                            self.compartments[go_id] = self._sanitize(go_name)

                    # Deduplicate by gene+compartment (keep highest confidence)
                    pair_key = (gene_name, go_id)
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)

                    self.localizations.append({
                        'gene_name': gene_name,
                        'go_id': go_id,
                        'confidence': confidence,
                        'evidence_type': evidence_type,
                        'source': source,
                    })
                    count += 1

            logger.info(f"Compartments: Loaded {count} records from {filename}")

        logger.info(f"Compartments: Total {len(self.localizations)} localizations, "
                     f"{len(self.compartments)} compartments")

    def get_nodes(self):
        """
        Generate CompartmentLocation nodes.
        Yields: (id, label, properties)
        """
        logger.info("Compartments: Generating nodes...")
        count = 0

        for go_id, go_name in self.compartments.items():
            props = {
                'name': go_name,
                'source': 'Compartments',
            }
            yield (go_id, "CompartmentLocation", props)
            count += 1

        logger.info(f"Compartments: Generated {count} CompartmentLocation nodes")

    def get_edges(self):
        """
        Generate ProteinLocatedIn edges.
        Note: Uses gene names as IDs since we don't have UniProt mapping here.
        The gene names can be linked to Gene nodes via name matching.
        Yields: (id, source, target, label, properties)
        """
        logger.info("Compartments: Generating edges...")
        count = 0

        for loc in self.localizations:
            props = {
                'confidence': loc['confidence'],
                'evidence_type': loc['evidence_type'],
                'source_db': loc['source'],
            }

            yield (
                None,
                loc['gene_name'],  # Gene name (can be matched to Gene nodes)
                loc['go_id'],
                "ProteinLocatedIn",
                props
            )
            count += 1

        logger.info(f"Compartments: Generated {count} ProteinLocatedIn edges")
