"""
OpenProt Adapter for BioCypher.

Loads OpenProt alternative protein data and generates:
- AlternativeProtein nodes (non-canonical ORFs and alt proteins)

OpenProt identifies novel proteins from alternative open reading frames,
including those from non-coding RNAs and alternative reading frames.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class OpenProtAdapter:
    def __init__(self, data_dir="template_package/data/openprot"):
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
        """Load OpenProt human summary data."""
        path = self.data_dir / 'openprot_human_summary.tsv'
        if not path.exists():
            logger.warning("OpenProt: summary file not found")
            return

        logger.info("OpenProt: Loading alternative protein data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                identifier = row.get('identifier', '').strip()
                prot_type = row.get('type', '').strip()
                gene_name = row.get('gene_name', '').strip()
                seq_length = row.get('sequence_length', '0').strip()
                transcript_acc = row.get('transcript_accessions', '').strip()
                protein_acc = row.get('protein_accessions', '').strip()

                if not identifier:
                    continue

                try:
                    length = int(seq_length)
                except ValueError:
                    length = 0

                self.proteins.append({
                    'id': identifier,
                    'type': prot_type,
                    'gene_name': gene_name,
                    'length': length,
                    'transcript_accessions': transcript_acc,
                    'protein_accessions': protein_acc,
                })
                count += 1

        logger.info(f"OpenProt: Loaded {count} protein entries")

    def get_nodes(self):
        """
        Generate AlternativeProtein nodes.
        Yields: (id, label, properties)
        """
        logger.info("OpenProt: Generating nodes...")
        count = 0

        for prot in self.proteins:
            props = {
                'protein_type': prot['type'],
                'gene_name': prot['gene_name'],
                'sequence_length': prot['length'],
                'transcript_accessions': self._sanitize(prot['transcript_accessions'][:200]),
                'source': 'OpenProt',
            }

            yield (prot['id'], "AlternativeProtein", props)
            count += 1

        logger.info(f"OpenProt: Generated {count} AlternativeProtein nodes")

    def get_edges(self):
        """No edges."""
        logger.info("OpenProt: No edges to generate")
        return iter([])
