"""
POSTAR3 Adapter for BioCypher.

Loads POSTAR3 post-transcriptional regulation data and generates:
- RBPClipIdentification edges (RBP genes identified by CLIP-seq methods)

POSTAR3 is a database for exploring post-transcriptional regulation.
The allTips.js file contains JavaScript arrays mapping RNA-binding protein
(RBP) gene names to the CLIP-seq methods that identified them (HITS-CLIP,
iCLIP, PAR-CLIP, eCLIP).

This adapter focuses on human RBP gene lists extracted from the JS data,
creating edges that link each gene symbol to the CLIP-seq method used to
identify it.
"""

import re
from pathlib import Path
from biocypher._logger import logger


# Mapping from JS variable name patterns to standardized CLIP method names.
# The allTips.js file uses variable names like human_clipdb_hits1,
# human_clipdb_i1, human_clipdb_par1, human_eclip to denote different
# CLIP-seq experimental methods.
CLIP_METHOD_MAP = {
    'hits': 'HITS-CLIP',
    'i': 'iCLIP',
    'par': 'PAR-CLIP',
    'eclip': 'eCLIP',
}


class POSTARAdapter:
    def __init__(self, data_dir="template_package/data/postar"):
        self.data_dir = Path(data_dir)
        self.clip_genes = {}  # {method_name: set(gene_symbols)}
        self._load_data()

    def _sanitize(self, text):
        """Clean and normalize a string value for safe output."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _classify_method(self, var_name):
        """
        Determine the CLIP-seq method from a JavaScript variable name.

        Variable naming conventions in allTips.js:
          human_clipdb_hits1, human_clipdb_hits2 -> HITS-CLIP
          human_clipdb_i1, human_clipdb_i2       -> iCLIP
          human_clipdb_par1, human_clipdb_par2   -> PAR-CLIP
          human_eclip                            -> eCLIP
          human_clipdb  (aggregate, no suffix)   -> skipped (composite list)

        Args:
            var_name: The JavaScript variable name (e.g. "human_clipdb_hits1").

        Returns:
            The standardized method name, or None if not classifiable.
        """
        if var_name == 'human_eclip':
            return 'eCLIP'

        # Match human_clipdb_{method_key}{number}
        match = re.match(r'^human_clipdb_([a-z]+)\d*$', var_name)
        if not match:
            return None

        key = match.group(1)
        return CLIP_METHOD_MAP.get(key)

    def _parse_js_arrays(self, content):
        """
        Parse JavaScript variable assignments containing string arrays.

        Each line in allTips.js has the form:
            variable_name=["val1","val2",...];

        Args:
            content: The full text content of the JS file.

        Returns:
            A list of (variable_name, list_of_strings) tuples.
        """
        results = []
        # Pattern: varname=["item1","item2",...];
        pattern = re.compile(
            r'^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*\[(.*?)\]\s*;?\s*$'
        )
        for line in content.splitlines():
            line = line.strip()
            if not line:
                continue
            match = pattern.match(line)
            if not match:
                continue
            var_name = match.group(1)
            array_body = match.group(2)
            # Extract quoted strings, handling escaped quotes
            items = re.findall(r'"([^"]*)"', array_body)
            results.append((var_name, items))
        return results

    def _load_data(self):
        """Load POSTAR3 RBP CLIP-seq gene lists from allTips.js."""
        path = self.data_dir / 'allTips.js'
        if not path.exists():
            logger.warning("POSTAR3: allTips.js not found at %s", path)
            return

        logger.info("POSTAR3: Loading CLIP-seq RBP data from allTips.js...")

        try:
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
        except Exception as e:
            logger.error("POSTAR3: Failed to read allTips.js: %s", e)
            return

        parsed = self._parse_js_arrays(content)

        for var_name, gene_list in parsed:
            method = self._classify_method(var_name)
            if method is None:
                continue

            if method not in self.clip_genes:
                self.clip_genes[method] = set()

            for gene in gene_list:
                gene = gene.strip()
                if gene:
                    self.clip_genes[method].add(gene)

        total = sum(len(genes) for genes in self.clip_genes.values())
        unique_genes = set()
        for genes in self.clip_genes.values():
            unique_genes.update(genes)
        logger.info(
            "POSTAR3: Loaded %d gene-method pairs across %d methods "
            "(%d unique gene symbols)",
            total,
            len(self.clip_genes),
            len(unique_genes),
        )

    def get_nodes(self):
        """
        Yield nodes. POSTAR3 RBPs reference existing Gene nodes, so no
        new nodes are created by this adapter.

        Yields:
            Nothing (empty iterator).
        """
        logger.info("POSTAR3: No new nodes (RBPs reference existing Gene nodes)")
        return iter([])

    def get_edges(self):
        """
        Yield RBPClipIdentification edges linking each RBP gene symbol
        to the CLIP-seq method that identified it.

        Each edge represents the fact that a given gene was identified as
        an RNA-binding protein by a specific CLIP-seq experimental method
        in the POSTAR3 database.

        Yields:
            Tuples of (id, source_id, target_id, label, properties_dict)
            where:
              - id: unique edge identifier (str)
              - source_id: gene symbol (str)
              - target_id: CLIP method identifier (str)
              - label: "RBPClipIdentification"
              - properties_dict: dict with source and method metadata
        """
        logger.info("POSTAR3: Generating RBPClipIdentification edges...")
        count = 0

        for method, genes in sorted(self.clip_genes.items()):
            method_id = f"CLIP:{method}"
            for gene in sorted(genes):
                edge_id = f"POSTAR3:{gene}:{method}"
                props = {
                    'gene_symbol': self._sanitize(gene),
                    'clip_method': method,
                    'source': 'POSTAR3',
                }
                yield (
                    edge_id,
                    gene,
                    method_id,
                    "RBPClipIdentification",
                    props,
                )
                count += 1

        logger.info("POSTAR3: Generated %d RBPClipIdentification edges", count)
