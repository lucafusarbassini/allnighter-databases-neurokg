"""
Cyclebase (Cell Cycle Regulation) Adapter for BioCypher.

Parses UniProt-derived TSV files containing cell cycle, cell division,
and CDK/cyclin protein annotations.

Yields CellCycleRegulation edges linking proteins to cell cycle phases/functions.
"""

from pathlib import Path
from biocypher._logger import logger


class CyclebaseAdapter:
    def __init__(self, data_dir="template_package/data/cyclebase"):
        self.data_dir = Path(data_dir)
        self.entries = []
        self._seen_keys = set()
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _extract_phase_from_keywords(self, keywords_str):
        """Extract cell cycle phase keywords from the Keywords field."""
        if not keywords_str:
            return ""
        kw_lower = keywords_str.lower()
        phases = []
        phase_terms = [
            "cell cycle", "cell division", "mitosis", "meiosis",
            "g1/s", "g2/m", "s phase", "m phase",
            "apoptosis", "dna damage", "dna repair",
            "kinetochore", "centromere", "chromosome",
        ]
        for term in phase_terms:
            if term in kw_lower:
                phases.append(term)
        return "; ".join(phases) if phases else "cell cycle"

    def _extract_function_summary(self, function_text):
        """Extract first sentence of FUNCTION annotation."""
        if not function_text:
            return ""
        text = function_text
        if text.startswith("FUNCTION: "):
            text = text[len("FUNCTION: "):]
        # Take first sentence (up to first period followed by space or end)
        dot_pos = text.find('. ')
        if dot_pos > 0:
            text = text[:dot_pos + 1]
        elif len(text) > 200:
            text = text[:200] + "..."
        return text

    def _classify_dataset(self, filename):
        """Return a category label based on the source filename."""
        fn = filename.lower()
        if "cdk_cyclin" in fn:
            return "CDK/cyclin"
        elif "celldivision" in fn:
            return "cell_division"
        elif "reviewed" in fn:
            return "cell_cycle_reviewed"
        else:
            return "cell_cycle"

    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("Cyclebase: data directory not found")
            return

        file_map = {
            "uniprot_cellcycle_human.tsv": "cell_cycle",
            "uniprot_cellcycle_reviewed.tsv": "cell_cycle_reviewed",
            "uniprot_celldivision_human.tsv": "cell_division",
            "uniprot_cdk_cyclins_human.tsv": "CDK/cyclin",
        }

        for fname, dataset in file_map.items():
            fpath = self.data_dir / fname
            if not fpath.exists():
                logger.warning(f"Cyclebase: {fname} not found, skipping")
                continue
            try:
                self._parse_tsv(fpath, dataset)
            except Exception as e:
                logger.warning(f"Cyclebase: Error reading {fpath}: {e}")

        logger.info(f"Cyclebase: Loaded {len(self.entries)} unique entries")

    def _parse_tsv(self, fpath, dataset):
        """Parse a UniProt-format TSV file."""
        with open(fpath, 'r', errors='replace') as fh:
            header_line = fh.readline().strip()
            if not header_line or header_line.startswith('<'):
                return
            headers = header_line.split('\t')

            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                row = dict(zip(headers, parts))

                uniprot_id = row.get('Entry', '').strip()
                if not uniprot_id:
                    continue

                # Extract gene name (first gene in the "Gene Names" field)
                gene_names_raw = row.get('Gene Names', '').strip()
                gene_name = gene_names_raw.split()[0] if gene_names_raw else ""

                # Dedup by (uniprot_id, dataset)
                dedup_key = (uniprot_id, dataset)
                if dedup_key in self._seen_keys:
                    continue
                self._seen_keys.add(dedup_key)

                keywords = row.get('Keywords', '').strip()
                phase = self._extract_phase_from_keywords(keywords)
                go_bp = row.get('Gene Ontology (biological process)', '').strip()
                function_cc = row.get('Function [CC]', '').strip()
                func_summary = self._extract_function_summary(function_cc)
                protein_name = row.get('Protein names', '').strip()

                self.entries.append({
                    'uniprot_id': uniprot_id,
                    'gene_name': gene_name,
                    'protein_name': protein_name,
                    'phase': phase,
                    'function': func_summary,
                    'keywords': keywords,
                    'go_bp': go_bp,
                    'dataset': dataset,
                })

    def get_nodes(self):
        logger.info("Cyclebase: No dedicated nodes (proteins referenced by UniProt ID)")
        return
        yield

    def get_edges(self):
        logger.info("Cyclebase: Generating CellCycleRegulation edges...")
        count = 0
        for entry in self.entries:
            edge_id = f"cyclebase:{entry['uniprot_id']}_{entry['dataset']}"
            source_id = entry['uniprot_id']
            target_id = f"cellcycle:{entry['dataset']}"

            props = {
                'phase': self._sanitize(entry['phase']),
                'function': self._sanitize(entry['function']),
                'gene_name': self._sanitize(entry['gene_name']),
                'protein_name': self._sanitize(entry['protein_name']),
                'keywords': self._sanitize(entry['keywords']),
                'go_biological_process': self._sanitize(entry['go_bp']),
                'dataset': entry['dataset'],
                'source': 'Cyclebase',
            }
            yield (edge_id, source_id, target_id,
                   "CellCycleRegulation", props)
            count += 1
        logger.info(f"Cyclebase: Generated {count} CellCycleRegulation edges")
