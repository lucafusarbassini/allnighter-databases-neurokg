"""
MitoCarta Adapter for BioCypher.

Loads MitoCarta 3.0 mitochondrial proteome data and generates:
- MitochondrialProtein edges (proteins localized to mitochondria)

MitoCarta is a comprehensive inventory of mitochondrial proteins,
scored by evidence from proteomics, literature, and homology.
"""

from html.parser import HTMLParser
from pathlib import Path
from biocypher._logger import logger


class _TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_row = False
        self.in_cell = False
        self.current_row = []
        self.rows = []

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.in_table = True
        elif tag == 'tr':
            self.in_row = True
            self.current_row = []
        elif tag in ('td', 'th'):
            self.in_cell = True

    def handle_endtag(self, tag):
        if tag == 'table':
            self.in_table = False
        elif tag == 'tr':
            self.in_row = False
            if self.current_row:
                self.rows.append(self.current_row)
        elif tag in ('td', 'th'):
            self.in_cell = False

    def handle_data(self, data):
        if self.in_cell:
            self.current_row.append(data.strip())


class MitoCartaAdapter:
    def __init__(self, data_dir="template_package/data/mitocarta"):
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
        """Load MitoCarta protein data from HTML."""
        path = self.data_dir / 'mitocarta_html.html'
        if not path.exists():
            logger.warning("MitoCarta: HTML data not found")
            return

        logger.info("MitoCarta: Loading mitochondrial proteome...")

        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()

        parser = _TableParser()
        parser.feed(content)

        if len(parser.rows) < 2:
            logger.warning("MitoCarta: No data rows found in HTML")
            return

        header = parser.rows[0]
        for row in parser.rows[1:]:
            if len(row) < 3:
                continue

            symbol = row[0].strip() if len(row) > 0 else ''
            description = row[1].strip() if len(row) > 1 else ''
            synonyms = row[2].strip() if len(row) > 2 else ''
            score = row[3].strip() if len(row) > 3 else '0'
            evidence = row[4].strip() if len(row) > 4 else ''

            if not symbol:
                continue

            # Extract UniProt ID from synonyms if present
            uniprot = ''
            for syn in synonyms.split(','):
                syn = syn.strip()
                if len(syn) == 6 and syn[0].isalpha() and syn[1:].replace('_', '').isalnum():
                    uniprot = syn
                    break

            try:
                score_int = int(score)
            except ValueError:
                score_int = 0

            self.proteins.append({
                'symbol': symbol,
                'description': description,
                'score': score_int,
                'evidence': evidence,
                'uniprot': uniprot,
            })

        logger.info(f"MitoCarta: Loaded {len(self.proteins)} mitochondrial proteins")

    def get_nodes(self):
        """No new nodes."""
        logger.info("MitoCarta: No new nodes")
        return iter([])

    def get_edges(self):
        """
        Generate MitochondrialProtein edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("MitoCarta: Generating edges...")
        count = 0

        for prot in self.proteins:
            props = {
                'gene_symbol': self._sanitize(prot['symbol']),
                'description': self._sanitize(prot['description'][:200]),
                'maestro_score': prot['score'],
                'evidence': self._sanitize(prot['evidence'][:200]),
                'source': 'MitoCarta3',
            }

            yield (
                None,
                prot['symbol'],
                "MITOCHONDRIA",
                "MitochondrialProtein",
                props
            )
            count += 1

        logger.info(f"MitoCarta: Generated {count} MitochondrialProtein edges")
