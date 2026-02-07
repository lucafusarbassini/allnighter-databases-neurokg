"""
CORUM Adapter for BioCypher.

Loads manually curated mammalian protein complex data from the CORUM database
(Comprehensive Resource of Mammalian Protein Complexes) and generates:
- ProteinComplex nodes (reuses the ComplexPortal schema type)
- ComplexContainsProtein edges (complex -> protein subunit membership)

Expected data directory: template_package/data/corum/
Expected file: allComplexes.txt or coreComplexes.txt (TSV with header)

CORUM TSV columns typically include:
ComplexID, ComplexName, Organism, Synonyms, Cell line, subunits(UniProt IDs),
subunits(Entrez IDs), Protein complex purification method, GO ID,
GO description, FunCat ID, FunCat description, subunits(Gene name),
Disease comment, PMID
"""

import zipfile
import io
from pathlib import Path
from biocypher._logger import logger


class CORUMAdapter:
    def __init__(self, data_dir="template_package/data/corum"):
        self.data_dir = Path(data_dir)
        self.complexes = {}
        self.subunits = []
        self._load_data()

    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _is_human(self, organism):
        """Check if the organism string refers to human."""
        if not organism:
            return False
        org_lower = organism.lower()
        return 'human' in org_lower or 'homo sapiens' in org_lower

    def _parse_subunits(self, subunits_str):
        """
        Parse semicolon-delimited UniProt subunit IDs.
        Example: 'P84022;Q13485;Q15796' -> ['P84022', 'Q13485', 'Q15796']
        Handles None, empty strings, and whitespace.
        """
        if not subunits_str or subunits_str.strip() in ('', '-', 'None'):
            return []
        parts = []
        for part in subunits_str.split(';'):
            uid = part.strip()
            if uid and uid != '-' and uid != 'None':
                parts.append(uid)
        return parts

    def _load_data(self):
        """Load CORUM complex data from TSV files, including inside zip archives."""
        if not self.data_dir.exists():
            logger.warning(f"CORUM: Data directory not found: {self.data_dir}")
            return

        files_parsed = 0

        # First try zip archives (CORUM sometimes distributes as .zip)
        for zf in self.data_dir.glob("*.zip"):
            try:
                with zipfile.ZipFile(zf, 'r') as z:
                    for name in z.namelist():
                        if name.endswith('.txt') or name.endswith('.tsv'):
                            with z.open(name) as f:
                                reader = io.TextIOWrapper(
                                    f, encoding='utf-8', errors='replace')
                                self._parse_tsv(reader)
                                files_parsed += 1
            except Exception as e:
                logger.warning(f"CORUM: Error reading {zf}: {e}")

        # Then try plain TSV/TXT files
        for fpath in (list(self.data_dir.glob("*.txt"))
                      + list(self.data_dir.glob("*.tsv"))):
            try:
                with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
                    first_line = f.readline()
                    # Skip HTML error pages or very short/corrupt files
                    if first_line.strip().startswith('<') or \
                            first_line.strip().startswith('404') or \
                            len(first_line.strip()) < 10:
                        logger.warning(
                            f"CORUM: {fpath.name} appears invalid, skipping")
                        continue
                    f.seek(0)
                    self._parse_tsv(f)
                    files_parsed += 1
            except Exception as e:
                logger.warning(f"CORUM: Error reading {fpath}: {e}")

        if files_parsed == 0:
            logger.warning(
                f"CORUM: No valid data files found in {self.data_dir}. "
                "Expected allComplexes.txt or coreComplexes.txt"
            )

        logger.info(
            f"CORUM: Loaded {len(self.complexes)} human complexes, "
            f"{len(self.subunits)} subunit links"
        )

    def _parse_tsv(self, fh):
        """
        Parse a CORUM TSV file from a file handle.
        Filters for human complexes and extracts complex + subunit data.
        """
        header = None
        for line in fh:
            if isinstance(line, bytes):
                line = line.decode('utf-8', errors='replace')
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            # First non-empty line is the header
            if header is None:
                header = parts
                continue

            if len(parts) < 5:
                continue

            row = dict(zip(header, parts))

            # Filter for human complexes only
            organism = row.get('Organism', row.get('organism', ''))
            if not self._is_human(organism):
                continue

            complex_id = row.get(
                'ComplexID', row.get('complex_id', '')).strip()
            if not complex_id:
                continue

            corum_id = f"CORUM:{complex_id}"
            complex_name = row.get(
                'ComplexName', row.get('complex_name', ''))

            # Extract subunit UniProt IDs
            subunits_str = row.get(
                'subunits(UniProt IDs)',
                row.get('subunits_uniprot_id', ''))
            parsed_subunits = self._parse_subunits(subunits_str)

            # Store complex info (deduplicated by ID)
            if corum_id not in self.complexes:
                self.complexes[corum_id] = {
                    'name': complex_name,
                    'organism': organism,
                    'synonyms': row.get('Synonyms', row.get('synonyms', '')),
                    'cell_line': row.get(
                        'Cell line', row.get('cell_line', '')),
                    'purification_method': row.get(
                        'Protein complex purification method',
                        row.get('purification_method', '')),
                    'go_id': row.get('GO ID', row.get('go_id', '')),
                    'go_description': row.get(
                        'GO description', row.get('go_description', '')),
                    'funcat_id': row.get(
                        'FunCat ID', row.get('funcat_id', '')),
                    'funcat_description': row.get(
                        'FunCat description',
                        row.get('funcat_description', '')),
                    'gene_names': row.get(
                        'subunits(Gene name)',
                        row.get('subunits_gene_name', '')),
                    'disease_comment': row.get(
                        'Disease comment', row.get('disease_comment', '')),
                    'pmid': row.get('PMID', row.get('pmid', '')),
                }

            # Store subunit links
            for uid in parsed_subunits:
                self.subunits.append({
                    'complex_id': corum_id,
                    'protein_id': uid,
                })

    def get_nodes(self):
        """
        Generate ProteinComplex nodes from CORUM data.
        Yields: (id, label, properties)
        """
        logger.info(
            f"CORUM: Generating nodes from {len(self.complexes)} complexes..."
        )
        count = 0

        for corum_id, info in self.complexes.items():
            num_components = sum(
                1 for s in self.subunits if s['complex_id'] == corum_id)

            props = {
                'name': self._sanitize(info['name']),
                'aliases': self._sanitize(info.get('synonyms', '')),
                'description': self._sanitize(
                    info.get('go_description', '')),
                'complex_assembly': '',
                'species': 'Homo sapiens',
                'go_annotations': [],
                'num_components': num_components,
                'go_id': self._sanitize(info.get('go_id', '')),
                'funcat_id': self._sanitize(info.get('funcat_id', '')),
                'funcat_description': self._sanitize(
                    info.get('funcat_description', '')),
                'gene_names': self._sanitize(info.get('gene_names', '')),
                'disease_comment': self._sanitize(
                    info.get('disease_comment', '')),
                'pmid': self._sanitize(info.get('pmid', '')),
                'purification_method': self._sanitize(
                    info.get('purification_method', '')),
                'cell_line': self._sanitize(info.get('cell_line', '')),
                'source': 'CORUM',
            }

            yield (corum_id, "ProteinComplex", props)
            count += 1

        logger.info(f"CORUM: Generated {count} ProteinComplex nodes")

    def get_edges(self):
        """
        Generate ComplexContainsProtein edges linking each complex to its
        protein subunits (UniProt IDs).
        Yields: (id, source, target, label, properties)
        """
        logger.info("CORUM: Generating ComplexContainsProtein edges...")
        count = 0
        unique_proteins = set()

        for sub in self.subunits:
            unique_proteins.add(sub['protein_id'])

            props = {
                'stoichiometry': 1,
                'species': 'Homo sapiens',
                'source': 'CORUM',
            }

            yield (
                None,
                sub['complex_id'],
                sub['protein_id'],
                "ComplexContainsProtein",
                props,
            )
            count += 1

        logger.info(
            f"CORUM: Generated {count} ComplexContainsProtein edges "
            f"(referencing {len(unique_proteins)} unique proteins)"
        )
