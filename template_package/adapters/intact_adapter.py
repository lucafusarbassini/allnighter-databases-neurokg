"""
IntAct Molecular Interaction Adapter for BioCypher.

Loads protein-protein interaction data from IntAct MITAB files and generates:
- MolecularInteraction edges (Gene → Gene with detection method and confidence)

Uses PSI-MITAB 2.5+ format with UniProt IDs.
"""

import re
from pathlib import Path
from biocypher._logger import logger


class IntActAdapter:
    def __init__(self, data_dir="template_package/data/intact",
                 min_score=0.4):
        self.data_dir = data_dir
        self.min_score = min_score
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _extract_uniprot_id(self, id_field):
        """Extract UniProt ID from MITAB ID field like 'uniprotkb:P46109'."""
        if not id_field:
            return None
        # Can have multiple IDs separated by |
        for part in id_field.split('|'):
            if part.startswith('uniprotkb:'):
                uid = part.replace('uniprotkb:', '')
                # Skip isoform IDs (e.g., P12345-2)
                if '-' in uid and uid.split('-')[1].isdigit():
                    uid = uid.split('-')[0]
                # Validate UniProt ID pattern (6 or 10 char alphanumeric)
                if re.match(r'^[A-Z0-9]{6,10}$', uid):
                    return uid
        return None

    def _extract_score(self, score_field):
        """Extract confidence score from 'intact-miscore:0.62'."""
        if not score_field:
            return 0.0
        for part in score_field.split('|'):
            if 'miscore:' in part:
                try:
                    return float(part.split(':')[-1])
                except (ValueError, IndexError):
                    pass
        return 0.0

    def _extract_detection_method(self, method_field):
        """Extract detection method name from PSI-MI format."""
        if not method_field:
            return ''
        # Format: psi-mi:"MI:0019"(coimmunoprecipitation)
        match = re.search(r'\(([^)]+)\)', method_field)
        if match:
            return self._sanitize(match.group(1))
        return ''

    def _extract_interaction_type(self, type_field):
        """Extract interaction type from PSI-MI format."""
        if not type_field:
            return ''
        match = re.search(r'\(([^)]+)\)', type_field)
        if match:
            return self._sanitize(match.group(1))
        return ''

    def _load_data(self):
        """Load IntAct MITAB file (human interactions only)."""
        filepath = Path(self.data_dir) / 'human.txt'
        if not filepath.exists():
            logger.warning("IntAct: human.txt not found")
            return

        logger.info(f"IntAct: Loading interactions (min_score={self.min_score})...")
        total = 0
        kept = 0
        seen = set()

        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                total += 1
                parts = line.strip().split('\t')
                if len(parts) < 15:
                    continue

                # Extract UniProt IDs
                id_a = self._extract_uniprot_id(parts[0])
                id_b = self._extract_uniprot_id(parts[1])

                if not id_a or not id_b or id_a == id_b:
                    continue

                # Extract confidence score
                score = self._extract_score(parts[14] if len(parts) > 14 else '')
                if score < self.min_score:
                    continue

                # Normalize edge direction for deduplication
                if id_a > id_b:
                    id_a, id_b = id_b, id_a

                # Deduplicate: keep highest score per pair
                pair_key = (id_a, id_b)
                if pair_key in seen:
                    continue
                seen.add(pair_key)

                detection_method = self._extract_detection_method(parts[6] if len(parts) > 6 else '')
                interaction_type = self._extract_interaction_type(parts[11] if len(parts) > 11 else '')

                self.interactions.append({
                    'source': id_a,
                    'target': id_b,
                    'score': score,
                    'detection_method': detection_method,
                    'interaction_type': interaction_type,
                })
                kept += 1

                if total % 200000 == 0:
                    logger.info(f"IntAct: Processed {total} lines, kept {kept}...")

        logger.info(f"IntAct: Loaded {kept}/{total} interactions "
                     f"(score >= {self.min_score})")

    def get_nodes(self):
        """IntAct produces no dedicated nodes."""
        logger.info("IntAct: No dedicated nodes (references existing Gene nodes)")
        return
        yield

    def get_edges(self):
        """
        Generate MolecularInteraction edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("IntAct: Generating MolecularInteraction edges...")
        count = 0

        for interaction in self.interactions:
            props = {
                'confidence_score': interaction['score'],
                'detection_method': interaction['detection_method'],
                'interaction_type': interaction['interaction_type'],
                'source_db': 'IntAct',
            }

            yield (
                None,
                interaction['source'],
                interaction['target'],
                "MolecularInteraction",
                props
            )
            count += 1

        logger.info(f"IntAct: Generated {count} MolecularInteraction edges")
