"""
NLSdb / UniProt Signal Peptide Adapter for BioCypher.

Loads UniProt signal peptide and transit peptide annotations and generates:
- SignalPeptide edges (proteins with signal/transit peptide annotations)

Provides localization signal information for human proteins.
"""

from pathlib import Path
from biocypher._logger import logger


class NLSdbAdapter:
    def __init__(self, data_dir="template_package/data/nlsdb"):
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
        """Load UniProt signal peptide data."""
        path = self.data_dir / 'uniprot_signals.tsv'
        if not path.exists():
            logger.warning("NLSdb: signal peptide data not found")
            return

        logger.info("NLSdb: Loading signal peptide annotations...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 4:
                    continue

                entry = parts[0].strip()
                gene_names = parts[1].strip()
                signal_peptide = parts[2].strip() if len(parts) > 2 else ''
                transit_peptide = parts[3].strip() if len(parts) > 3 else ''

                if not entry:
                    continue
                if not signal_peptide and not transit_peptide:
                    continue

                peptide_type = 'signal' if signal_peptide else 'transit'
                peptide_info = signal_peptide if signal_peptide else transit_peptide

                self.proteins.append({
                    'uniprot': entry,
                    'gene_name': gene_names.split()[0] if gene_names else '',
                    'peptide_type': peptide_type,
                    'peptide_info': peptide_info,
                })
                count += 1

        logger.info(f"NLSdb: Loaded {count} signal/transit peptide annotations")

    def get_nodes(self):
        """No new nodes."""
        logger.info("NLSdb: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate SignalPeptide edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("NLSdb: Generating edges...")
        seen = set()
        count = 0

        for prot in self.proteins:
            if prot['uniprot'] in seen:
                continue
            seen.add(prot['uniprot'])

            props = {
                'gene_name': self._sanitize(prot['gene_name']),
                'peptide_type': prot['peptide_type'],
                'peptide_info': self._sanitize(prot['peptide_info'][:200]),
                'source': 'UniProt_SignalP',
            }

            yield (
                None,
                prot['uniprot'],
                prot['peptide_type'],
                "SignalPeptide",
                props
            )
            count += 1

        logger.info(f"NLSdb: Generated {count} SignalPeptide edges")
