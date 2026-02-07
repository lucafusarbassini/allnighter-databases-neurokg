"""
PsychENCODE Adapter for BioCypher.

Loads PsychENCODE brain regulatory data and generates:
- EnhancerPromoterLink edges (enhancer → gene promoter Hi-C linkages)

PsychENCODE provides brain-specific regulatory element data including
Hi-C enhancer-promoter linkages and gene regulatory networks.
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class PsychENCODEAdapter:
    def __init__(self, data_dir="template_package/data/psychencode"):
        self.data_dir = Path(data_dir)
        self.links = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load PsychENCODE Hi-C enhancer-promoter linkages."""
        path = self.data_dir / 'INT-16_HiC_EP_linkages.csv'
        if not path.exists():
            logger.warning("PsychENCODE: Hi-C linkages file not found")
            return

        logger.info("PsychENCODE: Loading Hi-C enhancer-promoter linkages...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                chrom = row.get('Chromosome', '').strip()
                tss = row.get('Transcription_Start_Site', '').strip()
                gene_name = row.get('Target_Gene_Name', '').strip()
                ensembl = row.get('Target_Ensembl_Name', '').strip()
                enh_start = row.get('Enhancer_Start', '').strip()
                enh_end = row.get('Enhancer_End', '').strip()

                if not gene_name or not chrom:
                    continue

                try:
                    tss_int = int(tss)
                    enh_start_int = int(enh_start)
                    enh_end_int = int(enh_end)
                except ValueError:
                    continue

                self.links.append({
                    'chromosome': chrom,
                    'tss': tss_int,
                    'gene_name': gene_name,
                    'ensembl': ensembl,
                    'enh_start': enh_start_int,
                    'enh_end': enh_end_int,
                })
                count += 1

        logger.info(f"PsychENCODE: Loaded {count} enhancer-promoter linkages")

    def get_nodes(self):
        """No new nodes."""
        logger.info("PsychENCODE: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate EnhancerPromoterLink edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("PsychENCODE: Generating edges...")
        count = 0

        for link in self.links:
            enhancer_id = f"{link['chromosome']}:{link['enh_start']}-{link['enh_end']}"

            props = {
                'chromosome': link['chromosome'],
                'tss': link['tss'],
                'enhancer_start': link['enh_start'],
                'enhancer_end': link['enh_end'],
                'ensembl_id': link['ensembl'],
                'source': 'PsychENCODE',
            }

            yield (
                None,
                enhancer_id,
                link['gene_name'],
                "EnhancerPromoterLink",
                props
            )
            count += 1

        logger.info(f"PsychENCODE: Generated {count} EnhancerPromoterLink edges")
