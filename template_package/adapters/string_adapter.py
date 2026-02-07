"""
STRING Protein-Protein Interaction Adapter for BioCypher.

Loads STRING interaction data and generates:
- ProteinInteraction edges (Gene → Gene with confidence scores)

STRING uses Ensembl protein IDs (ENSP) which are mapped to UniProt IDs.
Only high-confidence interactions (combined_score >= 700) are included.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class STRINGAdapter:
    def __init__(self, data_dir="template_package/data/string",
                 min_score=700):
        self.data_dir = data_dir
        self.min_score = min_score
        self.ensp_to_uniprot = {}
        self.ensp_to_name = {}
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
        """Load STRING data files."""
        self._load_aliases()
        self._load_protein_info()
        self._load_interactions()

    def _load_aliases(self):
        """Load ENSP → UniProt ID mappings from aliases file."""
        filepath = Path(self.data_dir) / '9606.protein.aliases.v12.0.txt.gz'
        if not filepath.exists():
            logger.warning("STRING: aliases file not found")
            return

        logger.info("STRING: Loading UniProt ID mappings from aliases...")
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    ensp_id = parts[0]
                    alias = parts[1]
                    source = parts[2]

                    # Prefer UniProt_AC mappings (reviewed/canonical)
                    if source == 'UniProt_AC':
                        # Take the first UniProt AC as the primary one
                        # UniProt reviewed IDs are typically 6 chars
                        if ensp_id not in self.ensp_to_uniprot:
                            self.ensp_to_uniprot[ensp_id] = alias
                        elif len(alias) == 6 and len(self.ensp_to_uniprot[ensp_id]) != 6:
                            # Prefer canonical 6-char UniProt IDs
                            self.ensp_to_uniprot[ensp_id] = alias

        logger.info(f"STRING: Mapped {len(self.ensp_to_uniprot)} ENSP IDs to UniProt")

    def _load_protein_info(self):
        """Load protein names from info file."""
        filepath = Path(self.data_dir) / '9606.protein.info.v12.0.txt.gz'
        if not filepath.exists():
            return

        logger.info("STRING: Loading protein info...")
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    self.ensp_to_name[parts[0]] = parts[1]

    def _load_interactions(self):
        """Load protein-protein interactions above the score threshold."""
        filepath = Path(self.data_dir) / '9606.protein.links.v12.0.txt.gz'
        if not filepath.exists():
            logger.warning("STRING: links file not found")
            return

        logger.info(f"STRING: Loading interactions (min_score={self.min_score})...")
        total = 0
        kept = 0

        with gzip.open(filepath, 'rt') as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue

                total += 1
                protein1 = parts[0]
                protein2 = parts[1]
                score = int(parts[2])

                if score >= self.min_score:
                    # Map ENSP to UniProt
                    up1 = self.ensp_to_uniprot.get(protein1)
                    up2 = self.ensp_to_uniprot.get(protein2)

                    if up1 and up2 and up1 != up2:
                        # Normalize edge direction (alphabetical) for deduplication
                        if up1 > up2:
                            up1, up2 = up2, up1

                        self.interactions.append({
                            'source': up1,
                            'target': up2,
                            'score': score,
                            'ensp1': protein1,
                            'ensp2': protein2,
                        })
                        kept += 1

                if total % 2000000 == 0:
                    logger.info(f"STRING: Processed {total} links, kept {kept}...")

        logger.info(f"STRING: Loaded {kept}/{total} interactions "
                     f"(score >= {self.min_score})")

    def get_nodes(self):
        """
        STRING adapter produces no nodes directly.
        Interactions reference Gene nodes created by other adapters.
        Yields nothing.
        """
        logger.info("STRING: No dedicated nodes (references existing Gene nodes)")
        return
        yield  # Make it a generator

    def get_edges(self):
        """
        Generate ProteinInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("STRING: Generating ProteinInteraction edges...")
        count = 0

        for interaction in self.interactions:
            props = {
                'combined_score': interaction['score'],
                'source_db': 'STRING',
            }

            yield (
                None,
                interaction['source'],
                interaction['target'],
                "ProteinInteraction",
                props
            )
            count += 1

        logger.info(f"STRING: Generated {count} ProteinInteraction edges")
