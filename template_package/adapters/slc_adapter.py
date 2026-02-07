"""
SLC (Solute Carrier) Tables Adapter for BioCypher.

Loads SLC transporter protein data and generates:
- SLCTransporter nodes (SLC family members with function)

The SLC superfamily comprises over 400 membrane transport proteins
organized into 66 families.
"""

import csv
import json
from pathlib import Path
from biocypher._logger import logger


class SLCAdapter:
    def __init__(self, data_dir="template_package/data/slc"):
        self.data_dir = Path(data_dir)
        self.transporters = []
        self.families = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load SLC transporter data."""
        # Load UniProt TSV
        tsv_path = self.data_dir / 'slc_human_uniprot.tsv'
        if tsv_path.exists():
            logger.info("SLC: Loading transporter data from UniProt TSV...")
            with open(tsv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    entry = row.get('Entry', '').strip()
                    gene_names = row.get('Gene Names', '').strip()
                    protein_name = row.get('Protein names', '').strip()
                    function = row.get('Function [CC]', '').strip()

                    if not entry:
                        continue

                    # Extract primary gene name
                    primary_gene = gene_names.split()[0] if gene_names else ''

                    # Extract SLC family from gene name
                    family = ''
                    if primary_gene.upper().startswith('SLC'):
                        parts = primary_gene.upper().replace('SLC', '').split('A')
                        if parts:
                            family = f"SLC{parts[0]}"

                    # Truncate long function descriptions
                    func_short = function[:500] if function else ''

                    self.transporters.append({
                        'entry': entry,
                        'gene_name': primary_gene,
                        'protein_name': protein_name[:200],
                        'family': family,
                        'function': func_short,
                    })

            logger.info(f"SLC: Loaded {len(self.transporters)} transporter entries")

        # Load family summary
        fam_path = self.data_dir / 'slc_family_summary.json'
        if fam_path.exists():
            with open(fam_path, 'r', encoding='utf-8') as f:
                self.families = json.load(f)
            logger.info(f"SLC: Loaded {len(self.families)} family summaries")

    def get_nodes(self):
        """
        Generate SLCTransporter nodes.
        Yields: (id, label, properties)
        """
        logger.info("SLC: Generating nodes...")
        count = 0

        for t in self.transporters:
            props = {
                'gene_name': t['gene_name'],
                'protein_name': self._sanitize(t['protein_name']),
                'family': t['family'],
                'function': self._sanitize(t['function']),
                'source': 'SLC_Tables',
            }

            yield (t['entry'], "SLCTransporter", props)
            count += 1

        logger.info(f"SLC: Generated {count} SLCTransporter nodes")

    def get_edges(self):
        """No edges for now."""
        logger.info("SLC: No edges to generate")
        return iter([])
