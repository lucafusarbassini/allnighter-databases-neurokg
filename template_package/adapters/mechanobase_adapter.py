"""
MechanoBase Adapter for BioCypher.

Loads mechanosensitive protein annotation data from UniProt exports and generates:
- MechanosensitiveProtein edges (protein -> mechanosensitive annotation)

Data files:
- mechano_uniprot_v2.tsv       (50 core mechanosensitive proteins)
  Columns: Entry, Gene Names, Protein names, Function [CC], Keywords
- mechano_uniprot_channels.tsv  (477 mechanosensitive channel proteins)
  Columns: Entry, Entry Name, Gene Names, Organism, Protein names,
           Function [CC], Keywords, Gene Ontology (molecular function)

Proteins are annotated with mechanosensitive function from UniProt,
linking UniProt accession IDs to their functional annotations.
"""

from pathlib import Path
from biocypher._logger import logger


class MechanoBaseAdapter:
    def __init__(self, data_dir="template_package/data/mechanobase"):
        self.data_dir = Path(data_dir)
        self.core_proteins = []
        self.channel_proteins = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _truncate_function(self, text, max_length=500):
        """Truncate long function descriptions to avoid oversized properties."""
        text = self._sanitize(text)
        if len(text) > max_length:
            return text[:max_length] + '...'
        return text

    def _load_data(self):
        """Load both mechanosensitive protein data files."""
        self._load_core_proteins()
        self._load_channel_proteins()

    def _load_core_proteins(self):
        """Load mechano_uniprot_v2.tsv - core mechanosensitive proteins."""
        path = self.data_dir / 'mechano_uniprot_v2.tsv'
        if not path.exists():
            logger.warning("MechanoBase: mechano_uniprot_v2.tsv not found")
            return

        logger.info("MechanoBase: Loading core mechanosensitive proteins...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                entry = parts[0].strip()
                gene_names = parts[1].strip() if len(parts) > 1 else ''
                protein_names = parts[2].strip() if len(parts) > 2 else ''
                function_cc = parts[3].strip() if len(parts) > 3 else ''
                keywords = parts[4].strip() if len(parts) > 4 else ''

                if not entry:
                    continue

                # Extract primary gene name (first in space-separated list)
                primary_gene = gene_names.split()[0] if gene_names else ''

                self.core_proteins.append({
                    'uniprot_id': entry,
                    'gene_names': gene_names,
                    'primary_gene': primary_gene,
                    'protein_names': protein_names,
                    'function': function_cc,
                    'keywords': keywords,
                    'category': 'core_mechanosensitive',
                })
                count += 1

        logger.info(f"MechanoBase: Loaded {count} core mechanosensitive proteins")

    def _load_channel_proteins(self):
        """Load mechano_uniprot_channels.tsv - mechanosensitive channel proteins."""
        path = self.data_dir / 'mechano_uniprot_channels.tsv'
        if not path.exists():
            logger.warning("MechanoBase: mechano_uniprot_channels.tsv not found")
            return

        logger.info("MechanoBase: Loading mechanosensitive channel proteins...")
        count = 0

        with open(path, 'r', encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                entry = parts[0].strip()
                entry_name = parts[1].strip() if len(parts) > 1 else ''
                gene_names = parts[2].strip() if len(parts) > 2 else ''
                organism = parts[3].strip() if len(parts) > 3 else ''
                protein_names = parts[4].strip() if len(parts) > 4 else ''
                function_cc = parts[5].strip() if len(parts) > 5 else ''
                keywords = parts[6].strip() if len(parts) > 6 else ''
                go_mf = parts[7].strip() if len(parts) > 7 else ''

                if not entry:
                    continue

                primary_gene = gene_names.split()[0] if gene_names else ''

                self.channel_proteins.append({
                    'uniprot_id': entry,
                    'entry_name': entry_name,
                    'gene_names': gene_names,
                    'primary_gene': primary_gene,
                    'organism': organism,
                    'protein_names': protein_names,
                    'function': function_cc,
                    'keywords': keywords,
                    'go_molecular_function': go_mf,
                    'category': 'mechanosensitive_channel',
                })
                count += 1

        logger.info(f"MechanoBase: Loaded {count} mechanosensitive channel proteins")

    def get_nodes(self):
        """No new nodes - proteins are referenced by UniProt ID."""
        logger.info("MechanoBase: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate MechanosensitiveProtein annotation edges.

        Each edge links a UniProt protein ID to its mechanosensitive annotation.
        Source: UniProt ID, Target: primary gene name.

        Yields: (id, source, target, label, properties)
        """
        logger.info("MechanoBase: Generating edges...")
        count = 0
        seen = set()

        # Core mechanosensitive proteins
        for prot in self.core_proteins:
            if prot['uniprot_id'] in seen:
                continue
            seen.add(prot['uniprot_id'])

            props = {
                'gene_names': self._sanitize(prot['gene_names']),
                'protein_names': self._sanitize(prot['protein_names']),
                'function': self._truncate_function(prot['function']),
                'keywords': self._sanitize(prot['keywords']),
                'category': prot['category'],
                'source': 'MechanoBase',
            }

            yield (
                None,
                prot['uniprot_id'],
                prot['primary_gene'],
                "MechanosensitiveProtein",
                props
            )
            count += 1

        # Channel proteins
        for prot in self.channel_proteins:
            if prot['uniprot_id'] in seen:
                continue
            seen.add(prot['uniprot_id'])

            props = {
                'entry_name': prot['entry_name'],
                'gene_names': self._sanitize(prot['gene_names']),
                'organism': self._sanitize(prot['organism']),
                'protein_names': self._sanitize(prot['protein_names']),
                'function': self._truncate_function(prot['function']),
                'keywords': self._sanitize(prot['keywords']),
                'go_molecular_function': self._sanitize(prot['go_molecular_function']),
                'category': prot['category'],
                'source': 'MechanoBase',
            }

            yield (
                None,
                prot['uniprot_id'],
                prot['primary_gene'],
                "MechanosensitiveProtein",
                props
            )
            count += 1

        logger.info(f"MechanoBase: Generated {count} MechanosensitiveProtein edges")
