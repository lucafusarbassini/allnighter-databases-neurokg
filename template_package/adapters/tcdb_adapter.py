"""
TCDB (Transporter Classification Database) Adapter for BioCypher.

Loads TCDB transporter families and protein sequences, generating:
- TransporterFamily nodes (TCDB family classifications)
- ProteinInTransporterFamily edges (protein → family)

TCDB uses a hierarchical classification: class.subclass.family.subfamily.member
"""

import re
from pathlib import Path
from biocypher._logger import logger


class TCDBAdapter:
    def __init__(self, data_dir="template_package/data/tcdb"):
        self.data_dir = Path(data_dir)
        self.families = {}       # tc_id -> family_name
        self.proteins = []       # [{tc_id, pdb_or_uniprot, description}]
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load TCDB data files."""
        self._load_families()
        self._load_fasta()

    def _load_families(self):
        """Load TCDB family definitions from families.html (actually TSV)."""
        path = self.data_dir / 'families.html'
        if not path.exists():
            logger.warning("TCDB: families file not found")
            return

        logger.info("TCDB: Loading family definitions...")
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('<'):
                    continue
                parts = line.split('\t', 1)
                if len(parts) < 2:
                    continue
                tc_id = parts[0].strip()
                name = self._sanitize(parts[1])
                if tc_id and name:
                    self.families[tc_id] = name

        logger.info(f"TCDB: Loaded {len(self.families)} families")

    def _load_fasta(self):
        """Load TCDB FASTA file (tcdb.dat) for protein entries."""
        path = self.data_dir / 'tcdb.dat'
        if not path.exists():
            logger.warning("TCDB: FASTA file not found")
            return

        logger.info("TCDB: Loading protein entries from FASTA...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.startswith('>'):
                    continue

                # Parse FASTA header: >gnl|TC-DB|ACCESSION|TC_NUMBER Description
                header = line[1:].strip()
                parts = header.split('|')
                if len(parts) < 4:
                    continue

                accession = parts[2].strip()
                tc_and_desc = parts[3].strip()

                # Split TC number from description
                tc_parts = tc_and_desc.split(' ', 1)
                tc_number = tc_parts[0].strip()
                description = self._sanitize(tc_parts[1]) if len(tc_parts) > 1 else ''

                # Extract family-level TC ID (first 3 levels: X.Y.Z)
                tc_levels = tc_number.split('.')
                family_tc = '.'.join(tc_levels[:3]) if len(tc_levels) >= 3 else tc_number

                # Try to extract UniProt accession from the PDB-style accession
                # Some entries have UniProt IDs directly
                uniprot_id = ''
                if re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$', accession):
                    # Standard UniProt accession format
                    uniprot_id = accession
                elif re.match(r'^[A-Z0-9]+_[A-Z0-9]+$', accession):
                    # UniProt entry name (e.g., P0C1X8_HUMAN)
                    uniprot_id = accession

                self.proteins.append({
                    'accession': accession,
                    'tc_number': tc_number,
                    'family_tc': family_tc,
                    'description': description,
                    'uniprot_id': uniprot_id,
                })
                count += 1

        logger.info(f"TCDB: Loaded {count} protein entries")

    def get_nodes(self):
        """
        Generate TransporterFamily nodes.
        Yields: (id, label, properties)
        """
        logger.info("TCDB: Generating family nodes...")
        count = 0

        for tc_id, name in self.families.items():
            # Determine the class from the first digit
            tc_class = tc_id.split('.')[0] if '.' in tc_id else tc_id
            class_names = {
                '1': 'Channels/Pores',
                '2': 'Electrochemical Potential-driven Transporters',
                '3': 'Primary Active Transporters',
                '4': 'Group Translocators',
                '5': 'Transmembrane Electron Carriers',
                '8': 'Accessory Factors Involved in Transport',
                '9': 'Incompletely Characterized Transport Systems',
            }
            class_name = class_names.get(tc_class, '')

            props = {
                'name': name,
                'tc_class': class_name,
                'source': 'TCDB',
            }

            yield (f"TCDB:{tc_id}", "TransporterFamily", props)
            count += 1

        logger.info(f"TCDB: Generated {count} TransporterFamily nodes")

    def get_edges(self):
        """
        Generate ProteinInTransporterFamily edges.
        Links proteins (via UniProt or accession) to their TCDB family.
        Yields: (id, source, target, label, properties)
        """
        logger.info("TCDB: Generating edges...")
        seen = set()
        count = 0

        for prot in self.proteins:
            # Use UniProt ID if available, otherwise skip (no way to link to gene)
            source_id = prot['uniprot_id']
            if not source_id:
                continue

            family_tc = prot['family_tc']
            key = (source_id, family_tc)
            if key in seen:
                continue
            seen.add(key)

            # Only link if we have the family
            if family_tc not in self.families:
                continue

            props = {
                'tc_number': prot['tc_number'],
                'description': prot['description'],
                'source': 'TCDB',
            }

            yield (
                None,
                source_id,
                f"TCDB:{family_tc}",
                "ProteinInTransporterFamily",
                props
            )
            count += 1

        logger.info(f"TCDB: Generated {count} ProteinInTransporterFamily edges")
