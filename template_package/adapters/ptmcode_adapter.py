"""
PTMcode Adapter for BioCypher.

Loads PTMcode2 cross-protein PTM associations and generates:
- PTMCrosstalk edges (functional associations between PTMs on different proteins)

PTMcode identifies functional associations between post-translational
modifications across protein pairs using co-evolution, manual curation,
and 3D structural proximity.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class PTMcodeAdapter:
    def __init__(self, data_dir="template_package/data/ptmcode"):
        self.data_dir = Path(data_dir)
        self.associations = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PTMcode2 between-protein associations."""
        path = self.data_dir / 'PTMcode2_associations_between_proteins.txt.gz'
        if not path.exists():
            logger.warning("PTMcode: associations file not found")
            return

        logger.info("PTMcode: Loading cross-protein PTM associations...")
        count = 0
        skipped = 0

        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 13:
                    skipped += 1
                    continue

                prot1 = parts[0].strip()
                prot2 = parts[1].strip()
                species = parts[2].strip()

                # Filter for human
                if 'sapiens' not in species.lower() and 'human' not in species.lower():
                    # Also check for specific organism IDs
                    if species != 'Homo sapiens' and species != '9606':
                        skipped += 1
                        continue

                ptm1_type = parts[3].strip()
                residue1 = parts[4].strip()
                ptm2_type = parts[7].strip()
                residue2 = parts[8].strip()

                coevolution = parts[10].strip()
                manual = parts[11].strip()
                structure = parts[12].strip()

                self.associations.append({
                    'prot1': prot1,
                    'prot2': prot2,
                    'ptm1_type': ptm1_type,
                    'residue1': residue1,
                    'ptm2_type': ptm2_type,
                    'residue2': residue2,
                    'coevolution': coevolution == '1',
                    'manual': manual == '1',
                    'structure': structure == '1',
                })
                count += 1

                # Cap at 1M for memory
                if count >= 1000000:
                    break

        logger.info(f"PTMcode: Loaded {count} human cross-protein PTM associations (skipped {skipped})")

    def get_nodes(self):
        """No new nodes."""
        logger.info("PTMcode: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate PTMCrosstalk edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PTMcode: Generating edges...")
        seen = set()
        count = 0

        for assoc in self.associations:
            # Deduplicate by protein pair + residues
            key = (assoc['prot1'], assoc['prot2'], assoc['residue1'], assoc['residue2'])
            if key in seen:
                continue
            seen.add(key)

            evidence = []
            if assoc['coevolution']:
                evidence.append('coevolution')
            if assoc['manual']:
                evidence.append('manual')
            if assoc['structure']:
                evidence.append('structure')

            props = {
                'ptm1_type': assoc['ptm1_type'],
                'residue1': assoc['residue1'],
                'ptm2_type': assoc['ptm2_type'],
                'residue2': assoc['residue2'],
                'evidence_types': '|'.join(evidence),
                'source': 'PTMcode',
            }

            yield (
                None,
                assoc['prot1'],
                assoc['prot2'],
                "PTMCrosstalk",
                props
            )
            count += 1

        logger.info(f"PTMcode: Generated {count} PTMCrosstalk edges")
