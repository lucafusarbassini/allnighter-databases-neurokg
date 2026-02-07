"""
NuCOME (Nucleosome Organization Mapped across Epigenomes) Adapter for BioCypher.

Loads nucleosome positioning data (NPS = Nucleosome Positioning from Sequencing)
from CD34+ hematopoietic stem cells and generates:
- NucleosomePosition nodes (genomic positions of nucleosomes)

Data files:
- GSM651559_cd34_nucs_NPS.txt.gz  (2.5M nucleosome positions)

Columns: chr, start, end, summit, Most_likely_position, Tag_density,
         #Tag_positive_strand, #Tag_negative_strand, NPS_ID, -10(log10p-value)
"""

import gzip
from pathlib import Path
from biocypher._logger import logger


class NUCOMEAdapter:
    def __init__(self, data_dir="template_package/data/nucome"):
        self.data_dir = Path(data_dir)
        self.nucleosomes = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load nucleosome positioning data from NPS output."""
        path = self.data_dir / 'GSM651559_cd34_nucs_NPS.txt.gz'
        if not path.exists():
            logger.warning("NuCOME: GSM651559_cd34_nucs_NPS.txt.gz not found")
            return

        logger.info("NuCOME: Loading nucleosome positions from CD34+ cells...")
        count = 0

        with gzip.open(path, 'rt', encoding='utf-8') as f:
            # Skip header line
            header = f.readline()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 10:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                summit = int(parts[3])
                most_likely_pos = int(parts[4])
                tag_density = float(parts[5])
                tags_pos_strand = int(parts[6])
                tags_neg_strand = int(parts[7])
                nps_id = parts[8]
                neg10_log10_pval = float(parts[9])

                self.nucleosomes.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'summit': summit,
                    'most_likely_position': most_likely_pos,
                    'tag_density': tag_density,
                    'tags_positive_strand': tags_pos_strand,
                    'tags_negative_strand': tags_neg_strand,
                    'nps_id': nps_id,
                    'neg10_log10_pvalue': neg10_log10_pval,
                    'length': end - start,
                })
                count += 1

        logger.info(f"NuCOME: Loaded {count} nucleosome positions")

    def get_nodes(self):
        """
        Generate NucleosomePosition nodes.

        Yields: (id, label, properties)
          - id: NPS_ID (e.g., nucleosome960816)
        """
        logger.info("NuCOME: Generating NucleosomePosition nodes...")
        count = 0

        for nuc in self.nucleosomes:
            node_id = nuc['nps_id']
            props = {
                'chr': nuc['chrom'],
                'start': nuc['start'],
                'end': nuc['end'],
                'length': nuc['length'],
                'summit': nuc['summit'],
                'most_likely_position': nuc['most_likely_position'],
                'tag_density': nuc['tag_density'],
                'tags_positive_strand': nuc['tags_positive_strand'],
                'tags_negative_strand': nuc['tags_negative_strand'],
                'neg10_log10_pvalue': nuc['neg10_log10_pvalue'],
                'cell_type': 'CD34+',
                'source': 'NuCOME',
            }

            yield (node_id, "NucleosomePosition", props)
            count += 1

        logger.info(f"NuCOME: Generated {count} NucleosomePosition nodes")

    def get_edges(self):
        """No edges - nucleosome positions are standalone nodes."""
        logger.info("NuCOME: No edges")
        return iter([])
