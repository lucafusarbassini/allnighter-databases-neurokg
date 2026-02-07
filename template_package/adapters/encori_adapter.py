"""
ENCORI (starBase) Adapter for BioCypher.

Loads miRNA-target interaction data from ENCORI/starBase.
Data includes experimentally validated and computationally predicted
miRNA-mRNA interactions from CLIP-Seq experiments.

Generates:
- miRNA nodes (miRNA entities)
- miRNA-target edges (miRNA targeting gene relationships)
"""

import csv
from pathlib import Path
from biocypher._logger import logger


class ENCORIAdapter:
    def __init__(self, data_dir="template_package/data/encori"):
        self.data_dir = Path(data_dir)
        self.mirna_targets = []
        self.mirnas = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load ENCORI miRNA-target TSV files."""
        logger.info("ENCORI: Loading miRNA-target data...")
        total = 0

        for tsv_file in sorted(self.data_dir.glob('*.json')):
            # These are actually TSV files with .json extension and comment headers
            try:
                count = 0
                with open(tsv_file, 'r', encoding='utf-8') as f:
                    # Skip comment lines
                    header = None
                    for line in f:
                        line = line.strip()
                        if line.startswith('#'):
                            continue
                        if header is None:
                            header = line.split('\t')
                            continue

                        parts = line.split('\t')
                        if len(parts) < 10:
                            continue

                        row = dict(zip(header, parts))

                        mirna_id = row.get('miRNAid', '').strip()
                        mirna_name = row.get('miRNAname', '').strip()
                        gene_id = row.get('geneID', '').strip()
                        gene_name = row.get('geneName', '').strip()
                        gene_type = row.get('geneType', '').strip()
                        chrom = row.get('chromosome', '').strip()
                        clip_exp_num = row.get('clipExpNum', '0').strip()
                        target_scan = row.get('TargetScan', '0').strip()

                        if not mirna_id or not gene_id:
                            continue

                        # Only keep protein_coding gene targets
                        if gene_type and gene_type != 'protein_coding':
                            continue

                        self.mirnas[mirna_id] = mirna_name

                        clip_num = int(clip_exp_num) if clip_exp_num.isdigit() else 0
                        ts = int(target_scan) if target_scan.isdigit() else 0

                        self.mirna_targets.append({
                            'mirna_id': mirna_id,
                            'mirna_name': mirna_name,
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'chromosome': chrom,
                            'clip_experiments': clip_num,
                            'targetscan': ts,
                        })
                        count += 1

                logger.info(f"ENCORI: Loaded {count} targets from {tsv_file.name}")
                total += count
            except Exception as e:
                logger.warning(f"ENCORI: Error loading {tsv_file.name}: {e}")

        logger.info(f"ENCORI: Loaded {total} miRNA-target interactions total, {len(self.mirnas)} unique miRNAs")

    def get_nodes(self):
        """
        Generate miRNA nodes.
        Yields: (id, label, properties)
        """
        logger.info("ENCORI: Generating miRNA nodes...")
        count = 0

        for mirna_id, mirna_name in self.mirnas.items():
            props = {
                'name': self._sanitize(mirna_name),
                'source': 'ENCORI',
            }
            yield (mirna_id, "MicroRNA", props)
            count += 1

        logger.info(f"ENCORI: Generated {count} miRNA nodes")

    def get_edges(self):
        """
        Generate miRNA-target edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info(f"ENCORI: Generating edges from {len(self.mirna_targets)} interactions...")
        count = 0
        seen = set()

        for target in self.mirna_targets:
            pair = (target['mirna_id'], target['gene_id'])
            if pair in seen:
                continue
            seen.add(pair)

            props = {
                'gene_name': self._sanitize(target['gene_name']),
                'clip_experiments': target['clip_experiments'],
                'targetscan': target['targetscan'],
                'source': 'ENCORI',
            }

            yield (None, target['mirna_id'], target['gene_id'], "MiRNATargetsGene", props)
            count += 1

        logger.info(f"ENCORI: Generated {count} miRNA-target edges")
