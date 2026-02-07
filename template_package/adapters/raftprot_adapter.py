"""
RaftProt Adapter for BioCypher.

Loads lipid raft proteomics data from RaftProt (https://raftprot.org/).
RaftProt is a database of mammalian lipid raft-associated proteins
identified through detergent-resistant membrane (DRM) proteomics.

Generates:
- Gene-in-lipid-raft edges (proteins detected in lipid raft fractions)
"""

import csv
import json
from pathlib import Path
from biocypher._logger import logger


class RaftProtAdapter:
    def __init__(self, data_dir="template_package/data/raftprot"):
        self.data_dir = Path(data_dir)
        self.orthologs = {}
        self.proteins = []
        self._load_orthologs()
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_orthologs(self):
        """Load mouse-to-human ortholog mapping."""
        orth_path = Path("template_package/mappings/mouse_to_human_orthologs.json")
        if orth_path.exists():
            try:
                with open(orth_path, 'r') as f:
                    self.orthologs = json.load(f)
                logger.info(f"RaftProt: Loaded {len(self.orthologs)} mouse-to-human orthologs")
            except Exception as e:
                logger.warning(f"RaftProt: Could not load orthologs: {e}")

    def _load_data(self):
        """Load RaftProt data file."""
        data_path = self.data_dir / 'Raftprot.v2.4.txt'
        if not data_path.exists():
            logger.warning("RaftProt: Data file not found")
            return

        logger.info("RaftProt: Loading raft proteomics data...")
        count = 0

        try:
            with open(data_path, 'r', encoding='utf-8') as f:
                # File uses quoted, space-separated fields
                reader = csv.DictReader(f, delimiter=' ', quotechar='"')
                for row in reader:
                    uniprot = row.get('UniProt', row.get('OriginalID', '')).strip()
                    organism = row.get('Organism', '').strip()
                    gene_name = row.get('gene_name', '').strip()
                    protein_name = row.get('protein_name', '').strip()
                    method = row.get('BiochemMethod', '').strip()
                    detergent = row.get('Detergent', '').strip()
                    tissue_id = row.get('TissueID', '').strip()

                    if not uniprot:
                        continue

                    # Apply species filtering (Human + Mouse with ortholog mapping)
                    species = 'Homo sapiens' if organism == 'Human' else 'Mus musculus' if organism == 'Mouse' else organism
                    if species == 'Mus musculus':
                        mapped_id = self.orthologs.get(uniprot, uniprot)
                    elif species == 'Homo sapiens':
                        mapped_id = uniprot
                    else:
                        continue  # Skip other species

                    self.proteins.append({
                        'uniprot': mapped_id,
                        'original_id': uniprot,
                        'species': species,
                        'gene_name': gene_name,
                        'protein_name': protein_name,
                        'method': method,
                        'detergent': detergent,
                        'tissue_id': tissue_id,
                    })
                    count += 1
        except Exception as e:
            logger.warning(f"RaftProt: Error parsing data: {e}")

        logger.info(f"RaftProt: Loaded {count} raft protein entries")

    def get_nodes(self):
        """No additional nodes - uses existing Gene nodes."""
        logger.info("RaftProt: No additional nodes (uses Gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate gene-in-lipid-raft edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info(f"RaftProt: Generating edges from {len(self.proteins)} entries...")
        count = 0
        seen = set()

        # Create a single "lipid_raft" target node concept
        raft_target = "GO:0045121"  # lipid raft GO term

        for prot in self.proteins:
            uid = prot['uniprot']
            key = (uid, prot['species'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'gene_name': self._sanitize(prot['gene_name']),
                'method': self._sanitize(prot['method']),
                'detergent': self._sanitize(prot['detergent']),
                'species': prot['species'],
                'source': 'RaftProt',
            }

            if prot['species'] == 'Mus musculus' and prot['original_id'] != prot['uniprot']:
                props['original_id'] = prot['original_id']

            yield (None, uid, raft_target, "ProteinInLipidRaft", props)
            count += 1

        logger.info(f"RaftProt: Generated {count} raft association edges")
