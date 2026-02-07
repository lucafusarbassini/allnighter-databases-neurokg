"""
BioGRID (Biological General Repository for Interaction Datasets) Adapter for BioCypher.

Loads protein-protein interaction data from BioGRID tab3-format files and generates:
- ProteinInteraction edges (protein interacts with protein)

BioGRID is a curated database of protein and genetic interactions
from major model organisms. This adapter filters for human (taxid 9606)
physical and genetic interactions, using UniProt (SWISS-PROT) accessions
when available and falling back to Official Symbol identifiers.

Data format: BioGRID tab3 (TSV), typically distributed as a zip archive.
Expected file: BIOGRID-ALL-*.tab3.zip containing a single .txt TSV file.
"""

import zipfile
import io
from pathlib import Path
from biocypher._logger import logger

# BioGRID tab3 column indices (0-based)
COL_BIOGRID_ID = 0          # #BioGRID Interaction ID
COL_ENTREZ_A = 1             # Entrez Gene Interactor A
COL_ENTREZ_B = 2             # Entrez Gene Interactor B
COL_SYMBOL_A = 7             # Official Symbol Interactor A
COL_SYMBOL_B = 8             # Official Symbol Interactor B
COL_EXPERIMENTAL_SYSTEM = 11 # Experimental System
COL_EXPERIMENTAL_TYPE = 12   # Experimental System Type
COL_AUTHOR = 13              # Author
COL_PUBLICATION = 14         # Publication Source
COL_ORGANISM_A = 15          # Organism ID Interactor A
COL_ORGANISM_B = 16          # Organism ID Interactor B
COL_THROUGHPUT = 17          # Throughput
COL_SCORE = 18               # Score
COL_MODIFICATION = 19        # Modification
COL_SWISSPROT_A = 23         # SWISS-PROT Accessions Interactor A
COL_SWISSPROT_B = 26         # SWISS-PROT Accessions Interactor B

HUMAN_TAXID = '9606'


class BioGRIDAdapter:
    """
    Adapter for BioGRID protein-protein interaction data.

    Reads tab3-format TSV from zip archives, filters for human-human
    interactions, and yields edges with UniProt-based identifiers.
    """

    def __init__(self, data_dir="template_package/data/biogrid"):
        self.data_dir = Path(data_dir)
        self.interactions = []
        self._load_data()

    def _sanitize(self, text):
        """Clean text for safe CSV/property output."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _get_uniprot_id(self, swissprot_field, symbol):
        """
        Derive a protein identifier from available fields.

        Prefers SWISS-PROT accession (first accession if multiple are
        pipe-separated). Falls back to uniprot:{SYMBOL} format when
        no SWISS-PROT accession is available.

        Returns a string like 'P45985' or 'uniprot:MAP2K4'.
        """
        sp = swissprot_field.strip() if swissprot_field else ''
        if sp and sp != '-':
            # Take the first accession if multiple are listed (pipe-separated)
            first_acc = sp.split('|')[0].strip()
            if first_acc and first_acc != '-':
                return first_acc
        # Fall back to symbol-based identifier
        sym = symbol.strip() if symbol else ''
        if sym and sym != '-':
            return f"uniprot:{sym}"
        return None

    def _extract_pubmed_id(self, publication_field):
        """
        Extract PubMed ID from Publication Source field.

        Field is typically 'PUBMED:12345678'. Returns just the numeric
        ID, or the full field if format is unexpected.
        """
        pub = publication_field.strip() if publication_field else ''
        if pub.upper().startswith('PUBMED:'):
            return pub.split(':')[-1].strip()
        return pub

    def _find_zip_file(self):
        """Locate the BioGRID tab3 zip file in the data directory."""
        if not self.data_dir.exists():
            return None
        for entry in sorted(self.data_dir.iterdir()):
            if entry.suffix == '.zip' and 'tab3' in entry.name.lower():
                return entry
        # Fall back to any zip
        for entry in sorted(self.data_dir.iterdir()):
            if entry.suffix == '.zip' and 'biogrid' in entry.name.lower():
                return entry
        return None

    def _find_tab3_in_zip(self, zf):
        """Find the tab3 TSV file inside a zip archive."""
        for name in zf.namelist():
            if 'tab3' in name.lower() and name.endswith('.txt'):
                return name
        # Fall back to any .txt file
        txt_files = [n for n in zf.namelist() if n.endswith('.txt')]
        if txt_files:
            return txt_files[0]
        return None

    def _load_data(self):
        """
        Load BioGRID interactions from zip archive.

        Streams the file line by line to handle very large datasets
        (2.8M+ rows) without loading everything into memory at once.
        Filters for human-human interactions (both organisms = 9606).
        Deduplicates by (source, target, experimental_system, pubmed_id).
        """
        if not self.data_dir.exists():
            logger.warning("BioGRID: data directory not found")
            return

        zip_path = self._find_zip_file()
        if zip_path is None:
            logger.warning("BioGRID: no tab3 zip file found")
            return

        logger.info(f"BioGRID: Loading interactions from {zip_path.name}...")

        total_lines = 0
        human_count = 0
        kept = 0
        skipped_no_id = 0
        seen = set()

        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                tab3_name = self._find_tab3_in_zip(zf)
                if tab3_name is None:
                    logger.warning("BioGRID: no tab3 .txt file found inside zip")
                    return

                logger.info(f"BioGRID: Reading {tab3_name}...")

                with zf.open(tab3_name) as raw:
                    # Wrap in TextIOWrapper for line-by-line string reading
                    text_stream = io.TextIOWrapper(raw, encoding='utf-8', errors='replace')

                    # Skip header line
                    header = text_stream.readline()
                    if not header:
                        logger.warning("BioGRID: empty file")
                        return

                    for line in text_stream:
                        total_lines += 1

                        parts = line.split('\t')
                        if len(parts) < 27:
                            continue

                        # Filter for human-human interactions
                        org_a = parts[COL_ORGANISM_A].strip()
                        org_b = parts[COL_ORGANISM_B].strip()
                        if org_a != HUMAN_TAXID or org_b != HUMAN_TAXID:
                            continue

                        human_count += 1

                        # Extract identifiers
                        source_id = self._get_uniprot_id(
                            parts[COL_SWISSPROT_A], parts[COL_SYMBOL_A]
                        )
                        target_id = self._get_uniprot_id(
                            parts[COL_SWISSPROT_B], parts[COL_SYMBOL_B]
                        )

                        if not source_id or not target_id:
                            skipped_no_id += 1
                            continue

                        # Extract fields
                        biogrid_id = parts[COL_BIOGRID_ID].strip().lstrip('#')
                        experimental_system = parts[COL_EXPERIMENTAL_SYSTEM].strip()
                        experimental_type = parts[COL_EXPERIMENTAL_TYPE].strip()
                        throughput = parts[COL_THROUGHPUT].strip()
                        score_raw = parts[COL_SCORE].strip()
                        score = score_raw if score_raw != '-' else ''
                        modification = parts[COL_MODIFICATION].strip()
                        modification = modification if modification != '-' else ''
                        pubmed_id = self._extract_pubmed_id(parts[COL_PUBLICATION])

                        # Deduplicate by (source, target, exp_system, pubmed)
                        # Normalize edge direction for dedup
                        norm_src = min(source_id, target_id)
                        norm_tgt = max(source_id, target_id)
                        dedup_key = (norm_src, norm_tgt, experimental_system, pubmed_id)
                        if dedup_key in seen:
                            continue
                        seen.add(dedup_key)

                        self.interactions.append({
                            'biogrid_id': biogrid_id,
                            'source_id': source_id,
                            'target_id': target_id,
                            'symbol_a': parts[COL_SYMBOL_A].strip(),
                            'symbol_b': parts[COL_SYMBOL_B].strip(),
                            'experimental_system': experimental_system,
                            'experimental_type': experimental_type,
                            'throughput': throughput,
                            'score': score,
                            'modification': modification,
                            'pubmed_id': pubmed_id,
                        })
                        kept += 1

                        if total_lines % 500000 == 0:
                            logger.info(
                                f"BioGRID: Processed {total_lines} lines, "
                                f"{human_count} human, {kept} kept..."
                            )

        except (zipfile.BadZipFile, OSError) as exc:
            logger.error(f"BioGRID: Error reading zip file: {exc}")
            return

        logger.info(
            f"BioGRID: Loaded {kept} unique human interactions "
            f"from {human_count} human rows ({total_lines} total lines). "
            f"Skipped {skipped_no_id} rows with missing IDs."
        )

    def get_nodes(self):
        """
        No dedicated nodes -- proteins/genes are handled by the
        unified Gene entity from UniProt or other gene adapters.
        """
        logger.info("BioGRID: No dedicated nodes (references existing protein/gene nodes)")
        return
        yield

    def get_edges(self):
        """
        Generate protein-protein interaction edges.

        Yields tuples of:
            (edge_id, source_id, target_id, label, properties)

        Label: 'protein interacts with protein'
        """
        logger.info("BioGRID: Generating protein interaction edges...")
        count = 0

        for ix in self.interactions:
            props = {
                'biogrid_interaction_id': self._sanitize(ix['biogrid_id']),
                'symbol_a': self._sanitize(ix['symbol_a']),
                'symbol_b': self._sanitize(ix['symbol_b']),
                'experimental_system': self._sanitize(ix['experimental_system']),
                'experimental_type': self._sanitize(ix['experimental_type']),
                'throughput': self._sanitize(ix['throughput']),
                'pubmed_id': self._sanitize(ix['pubmed_id']),
                'score': self._sanitize(ix['score']),
                'modification': self._sanitize(ix['modification']),
                'source': 'BioGRID',
            }

            yield (
                None,
                ix['source_id'],
                ix['target_id'],
                'protein interacts with protein',
                props,
            )
            count += 1

        logger.info(f"BioGRID: Generated {count} 'protein interacts with protein' edges")
