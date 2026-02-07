"""
iPTMnet Adapter for BioCypher.

Loads iPTMnet post-translational modification data and generates:
- ProteinHasPTM edges (protein → PTM site with type, enzyme, and evidence)

iPTMnet integrates PTM information from multiple sources including
PRO, PhosphoSitePlus, and literature mining.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class iPTMnetAdapter:
    def __init__(self, data_dir="template_package/data/iptmnet"):
        self.data_dir = Path(data_dir)
        self.ptm_records = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load iPTMnet compact JSONL data."""
        path = self.data_dir / 'iptmnet_human_ptm_compact.jsonl'
        if not path.exists():
            logger.warning("iPTMnet: data file not found")
            return

        logger.info("iPTMnet: Loading PTM data...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:
                    continue

                accession = entry.get('accession', '')
                substrate_data = entry.get('substrate_data', {})

                if not accession or not substrate_data:
                    continue

                # Flatten all PTMs across isoforms
                for isoform_id, ptms in substrate_data.items():
                    for ptm in ptms:
                        site = ptm.get('site', '')
                        ptm_type = ptm.get('ptm_type', '')
                        residue = ptm.get('residue', '')
                        score = ptm.get('score', 0)

                        # Get enzyme info
                        enzymes = ptm.get('enzymes', [])
                        enzyme_names = [e.get('name', '') for e in enzymes if e.get('name')]

                        # Get PMIDs
                        pmids = ptm.get('pmids', [])

                        self.ptm_records.append({
                            'accession': accession,
                            'site': site,
                            'ptm_type': ptm_type,
                            'residue': residue,
                            'score': score,
                            'enzymes': '|'.join(enzyme_names),
                            'pmids': '|'.join(pmids),
                        })
                        count += 1

        logger.info(f"iPTMnet: Loaded {count} PTM records")

    def get_nodes(self):
        """No new nodes - links to existing Gene nodes via UniProt ID."""
        logger.info("iPTMnet: No new nodes (uses existing gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Generate ProteinHasPTM edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("iPTMnet: Generating edges...")
        seen = set()
        count = 0

        for rec in self.ptm_records:
            # Deduplicate by accession + site + type
            key = (rec['accession'], rec['site'], rec['ptm_type'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'site': rec['site'],
                'ptm_type': rec['ptm_type'],
                'residue': rec['residue'],
                'score': rec['score'],
                'enzymes': self._sanitize(rec['enzymes']),
                'pmids': rec['pmids'],
                'source': 'iPTMnet',
            }

            yield (
                None,
                rec['accession'],
                f"PTM:{rec['accession']}_{rec['site']}",
                "ProteinHasPTM",
                props
            )
            count += 1

        logger.info(f"iPTMnet: Generated {count} ProteinHasPTM edges")
