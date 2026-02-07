"""
BioLiP Adapter for BioCypher.

Loads BioLiP ligand-protein interaction data and generates:
- Ligand nodes (unique small molecule ligands from PDB)
- LigandBindsProtein edges (ligand → protein with binding site info)

Uses the non-redundant (nr) set to avoid duplicate structures.
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class BioLiPAdapter:
    def __init__(self, data_dir="template_package/data/biolip"):
        self.data_dir = Path(data_dir)
        self.ligands = {}          # ligand_id -> {name, ec_number, go_terms}
        self.interactions = []     # [{pdb_id, ligand_id, uniprot_id, ...}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load BioLiP non-redundant data."""
        nr_path = self.data_dir / 'BioLiP_nr.txt.gz'
        if not nr_path.exists():
            nr_path = self.data_dir / 'BioLiP_nr.txt'
            if not nr_path.exists():
                logger.warning("BioLiP: Data file not found")
                return

        logger.info("BioLiP: Loading non-redundant ligand-protein interactions...")
        count = 0
        skipped = 0

        open_func = gzip.open if str(nr_path).endswith('.gz') else open
        with open_func(nr_path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 18:
                    skipped += 1
                    continue

                pdb_id = parts[0].strip()
                receptor_chain = parts[1].strip()
                resolution = parts[2].strip()
                ligand_id = parts[4].strip()       # CCD ligand ID
                ligand_chain = parts[5].strip()
                ec_number = parts[11].strip()
                go_terms = parts[12].strip()
                uniprot_id = parts[17].strip() if len(parts) > 17 else ''

                # Skip DNA/RNA/peptide ligands (keep small molecules)
                if ligand_id.lower() in ('dna', 'rna', 'peptide', ''):
                    continue

                # Register ligand
                if ligand_id not in self.ligands:
                    self.ligands[ligand_id] = {
                        'ec_number': ec_number,
                        'go_terms': go_terms,
                    }

                # Register interaction (only if we have a UniProt ID)
                if uniprot_id and len(uniprot_id) >= 4:
                    self.interactions.append({
                        'pdb_id': pdb_id,
                        'receptor_chain': receptor_chain,
                        'resolution': resolution,
                        'ligand_id': ligand_id,
                        'ligand_chain': ligand_chain,
                        'uniprot_id': uniprot_id,
                        'ec_number': ec_number,
                    })
                    count += 1

        logger.info(f"BioLiP: Loaded {len(self.ligands)} unique ligands, "
                     f"{count} interactions (skipped {skipped} malformed lines)")

    def get_nodes(self):
        """
        Generate Ligand nodes (unique small molecule ligands).
        Yields: (id, label, properties)
        """
        logger.info("BioLiP: Generating ligand nodes...")
        count = 0

        for ligand_id, data in self.ligands.items():
            props = {
                'name': ligand_id,
                'ec_number': self._sanitize(data.get('ec_number', '')),
                'go_terms': self._sanitize(data.get('go_terms', '')),
                'source': 'BioLiP',
            }

            yield (f"PDB:{ligand_id}", "PDBLigand", props)
            count += 1

        logger.info(f"BioLiP: Generated {count} PDBLigand nodes")

    def get_edges(self):
        """
        Generate LigandBindsProtein edges.
        Deduplicated by (ligand, uniprot) pair.
        Yields: (id, source, target, label, properties)
        """
        logger.info("BioLiP: Generating edges...")
        seen = set()
        count = 0

        for inter in self.interactions:
            key = (inter['ligand_id'], inter['uniprot_id'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'pdb_id': inter['pdb_id'],
                'resolution': inter['resolution'],
                'ec_number': self._sanitize(inter['ec_number']),
                'source': 'BioLiP',
            }

            yield (
                None,
                f"PDB:{inter['ligand_id']}",
                inter['uniprot_id'],
                "LigandBindsProtein",
                props
            )
            count += 1

        logger.info(f"BioLiP: Generated {count} LigandBindsProtein edges")
