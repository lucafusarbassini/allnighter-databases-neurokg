"""
GlyGen (Glycan Database) Adapter for BioCypher.

Loads glycan structures from GlyGen/GlyTouCan and generates:
- Glycan nodes (glycan structures with mass and composition)

Data sourced from GlyGen API (api.glygen.org) which integrates
UniCarbKB, GlyTouCan, and other glycoinformatics resources.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class GlyGenAdapter:
    def __init__(self, data_dir="template_package/data/unicarbkb"):
        self.data_dir = Path(data_dir)
        self.glycans = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load GlyGen glycan TSV files."""
        # Load N-linked glycans (larger dataset)
        nlinked_path = self.data_dir / 'glygen_nlinked_glycans.tsv'
        if nlinked_path.exists():
            self._load_tsv(nlinked_path, 'N-linked')

        # Load general glycans
        general_path = self.data_dir / 'glygen_glycans.tsv'
        if general_path.exists():
            self._load_tsv(general_path, 'general')

        # Deduplicate by GlyTouCan accession
        seen = set()
        unique = []
        for g in self.glycans:
            if g['glytoucan_ac'] not in seen:
                seen.add(g['glytoucan_ac'])
                unique.append(g)
        self.glycans = unique

        logger.info(f"GlyGen: Loaded {len(self.glycans)} unique glycan structures")

    def _load_tsv(self, path, glycan_type):
        """Load a TSV file of glycan records."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    acc = row.get('glytoucan_ac', '').strip()
                    if not acc:
                        continue

                    mass = row.get('mass', '')
                    byonic = row.get('byonic', '')
                    # Parse composition from byonic field
                    composition = byonic.split('%')[0].strip() if byonic and '%' in byonic else byonic

                    self.glycans.append({
                        'glytoucan_ac': acc,
                        'mass': float(mass) if mass else 0.0,
                        'composition': self._sanitize(composition),
                        'glycan_type': glycan_type,
                        'hit_score': float(row.get('hit_score', 0)) if row.get('hit_score') else 0.0,
                        'publication_count': int(row.get('publication_count', 0)) if row.get('publication_count') else 0,
                    })
        except Exception as e:
            logger.warning(f"GlyGen: Error loading {path}: {e}")

    def get_nodes(self):
        """
        Generate Glycan nodes.
        Yields: (id, label, properties)
        """
        logger.info("GlyGen: Generating glycan nodes...")
        count = 0

        for glycan in self.glycans:
            props = {
                'composition': glycan['composition'],
                'mass': glycan['mass'],
                'glycan_type': glycan['glycan_type'],
                'hit_score': glycan['hit_score'],
                'publication_count': glycan['publication_count'],
                'source': 'GlyGen',
            }

            yield (f"GlyTouCan:{glycan['glytoucan_ac']}", "Glycan", props)
            count += 1

        logger.info(f"GlyGen: Generated {count} Glycan nodes")

    def get_edges(self):
        """
        No edges for glycan nodes alone.
        Glycan-protein links would require glycoprotein annotation data.
        """
        logger.info("GlyGen: No edges to generate")
        return iter([])
