"""
DrLLPS Adapter for BioCypher.

Loads DrLLPS liquid-liquid phase separation data and generates:
- LLPSProtein edges (proteins involved in LLPS with condensate info)

DrLLPS is a database of proteins driving liquid-liquid phase separation,
providing information on scaffolds, clients, and regulators of condensates.
"""

from pathlib import Path
from biocypher._logger import logger


class DrLLPSAdapter:
    def __init__(self, data_dir="template_package/data/drllps"):
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
        """Load DrLLPS LLPS protein data."""
        path = self.data_dir / 'drllps_all.txt'
        if not path.exists():
            logger.warning("DrLLPS: data file not found")
            return

        logger.info("DrLLPS: Loading LLPS protein data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = None
            for line in f:
                parts = line.strip().split('\t')
                if header is None:
                    header = parts
                    continue

                if len(parts) < 7:
                    continue

                drllps_id = parts[0].strip()
                uniprot = parts[1].strip()
                gene_name = parts[2].strip()
                ensembl = parts[3].strip()
                species = parts[4].strip()
                condensate = parts[5].strip()
                llps_type = parts[6].strip()

                if species != 'Homo sapiens':
                    continue

                if not uniprot:
                    continue

                self.proteins.append({
                    'drllps_id': drllps_id,
                    'uniprot': uniprot,
                    'gene_name': gene_name,
                    'ensembl': ensembl,
                    'condensate': condensate,
                    'llps_type': llps_type,
                })
                count += 1

        logger.info(f"DrLLPS: Loaded {count} human LLPS proteins")

    def get_nodes(self):
        """No new nodes."""
        logger.info("DrLLPS: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate LLPSProtein edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("DrLLPS: Generating edges...")
        seen = set()
        count = 0

        for prot in self.proteins:
            if prot['uniprot'] in seen:
                continue
            seen.add(prot['uniprot'])

            props = {
                'gene_name': self._sanitize(prot['gene_name']),
                'condensate': self._sanitize(prot['condensate']),
                'llps_type': prot['llps_type'],
                'ensembl_id': prot['ensembl'],
                'source': 'DrLLPS',
            }

            yield (
                None,
                prot['uniprot'],
                "LLPS",
                "LLPSProtein",
                props
            )
            count += 1

        logger.info(f"DrLLPS: Generated {count} LLPSProtein edges")
