"""
GWAS Catalog Adapter for BioCypher.
Loads genome-wide association study data from the EBI GWAS Catalog.

Generates:
- Disease/Trait nodes from EFO ontology terms in the mapped trait URIs
- GWASAssociation edges linking SNP variants to diseases/traits

Data files:
- gwas-catalog-associations_ontology-annotated-full.zip (associations TSV)
- gwas-catalog-studies.tsv (study metadata, keyed by PUBMEDID)

Only genome-wide significant associations (p-value < 5e-8) are included.
"""

import csv
import io
import zipfile
from pathlib import Path
from biocypher._logger import logger


class GWASCatalogAdapter:
    # Standard GWAS significance threshold
    PVALUE_THRESHOLD = 5e-8

    def __init__(self, data_dir="template_package/data/gwas_catalog"):
        self.data_dir = Path(data_dir)
        self.traits = {}         # EFO URI -> trait info
        self.associations = []   # significant variant-trait associations
        self.studies = {}        # pubmed_id -> study metadata
        self._load_studies()
        self._load_associations()

    def _sanitize(self, text):
        """Clean text for safe CSV/property storage."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _parse_pvalue(self, pval_str):
        """
        Parse p-value string to float.
        Returns None if unparseable.
        Handles scientific notation like '9E-24', '3E-19', etc.
        """
        if not pval_str or not pval_str.strip():
            return None
        try:
            return float(pval_str.strip())
        except (ValueError, TypeError):
            return None

    def _extract_efo_id(self, uri):
        """
        Extract an EFO-style ID from a trait URI.
        e.g., 'http://www.ebi.ac.uk/efo/EFO_0004713' -> 'EFO:0004713'
              'http://purl.obolibrary.org/obo/HP_0001250' -> 'HP:0001250'
        """
        if not uri:
            return None
        uri = uri.strip()
        # Extract the last path component
        fragment = uri.rstrip('/').rsplit('/', 1)[-1]
        # Convert underscore to colon for standard ontology ID format
        if '_' in fragment:
            parts = fragment.split('_', 1)
            return f"{parts[0]}:{parts[1]}"
        return fragment

    def _load_studies(self):
        """
        Load study metadata from the studies TSV file.
        The studies file is keyed by PUBMEDID (no STUDY ACCESSION column).
        We store one entry per PUBMEDID for enriching association edges.
        """
        path = self.data_dir / 'gwas-catalog-studies.tsv'
        if not path.exists():
            logger.warning(
                "GWAS Catalog: studies file not found, "
                "proceeding without study metadata"
            )
            return

        logger.info("GWAS Catalog: Loading study metadata...")
        count = 0

        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pubmed_id = (row.get('PUBMEDID') or '').strip()
                if not pubmed_id:
                    continue

                # Store first entry per PUBMEDID (some have multiple traits)
                if pubmed_id not in self.studies:
                    self.studies[pubmed_id] = {
                        'title': (row.get('STUDY') or '').strip(),
                        'disease_trait': (row.get('DISEASE/TRAIT') or '').strip(),
                        'initial_sample': (
                            row.get('INITIAL SAMPLE SIZE') or ''
                        ).strip(),
                        'first_author': (row.get('FIRST AUTHOR') or '').strip(),
                        'journal': (row.get('JOURNAL') or '').strip(),
                        'date': (row.get('DATE') or '').strip(),
                    }
                    count += 1

        logger.info(f"GWAS Catalog: Loaded metadata for {count} studies")

    def _load_associations(self):
        """
        Load GWAS associations from the zipped TSV.

        Filters for:
        - Valid SNP ID (rsXXXXXX format)
        - Valid mapped trait URI
        - Genome-wide significance (p-value < 5e-8)

        Deduplicates by (SNP, trait_id) pairs, keeping the most
        significant association (lowest p-value).
        """
        zip_path = (
            self.data_dir
            / 'gwas-catalog-associations_ontology-annotated-full.zip'
        )
        if not zip_path.exists():
            logger.warning("GWAS Catalog: associations zip file not found")
            return

        logger.info("GWAS Catalog: Loading associations (p < 5e-8)...")
        skipped_pval = 0
        skipped_snp = 0
        skipped_trait = 0

        # Track best association per (snp, trait_id) pair
        best = {}

        with zipfile.ZipFile(zip_path, 'r') as zf:
            # Find the TSV inside the zip
            tsv_names = [n for n in zf.namelist() if n.endswith('.tsv')]
            if not tsv_names:
                logger.warning("GWAS Catalog: no TSV found inside zip")
                return

            tsv_name = tsv_names[0]
            logger.info(f"GWAS Catalog: Reading {tsv_name} from zip...")

            with zf.open(tsv_name) as raw:
                text_stream = io.TextIOWrapper(
                    raw, encoding='utf-8', errors='replace'
                )
                reader = csv.DictReader(text_stream, delimiter='\t')

                for row in reader:
                    # Extract SNP ID
                    snps = (row.get('SNPS') or '').strip()
                    if not snps or not snps.startswith('rs'):
                        skipped_snp += 1
                        continue

                    # Handle multiple SNPs: take first rs ID
                    snp_id = snps.split(';')[0].strip().split(',')[0].strip()
                    if not snp_id.startswith('rs'):
                        for part in snps.replace(';', ' ').replace(',', ' ').split():
                            if part.startswith('rs'):
                                snp_id = part
                                break
                        else:
                            skipped_snp += 1
                            continue

                    # Parse p-value and filter for significance
                    pval = self._parse_pvalue(row.get('P-VALUE', ''))
                    if pval is None or pval >= self.PVALUE_THRESHOLD:
                        skipped_pval += 1
                        continue

                    # Get mapped trait URIs (can be comma-separated)
                    trait_uri_raw = (
                        row.get('MAPPED_TRAIT_URI') or ''
                    ).strip()
                    trait_name_raw = (row.get('MAPPED_TRAIT') or '').strip()
                    if not trait_uri_raw:
                        skipped_trait += 1
                        continue

                    # Split multiple traits
                    trait_uris = [
                        u.strip() for u in trait_uri_raw.split(',')
                        if u.strip()
                    ]
                    trait_names = [
                        n.strip() for n in trait_name_raw.split(',')
                        if n.strip()
                    ]

                    # Parse association properties
                    risk_allele = (
                        row.get('STRONGEST SNP-RISK ALLELE') or ''
                    ).strip()
                    or_beta_str = (row.get('OR or BETA') or '').strip()
                    ci_text = (row.get('95% CI (TEXT)') or '').strip()
                    pubmed_id = (row.get('PUBMEDID') or '').strip()
                    study_accession = (
                        row.get('STUDY ACCESSION') or ''
                    ).strip()
                    mapped_gene = (row.get('MAPPED_GENE') or '').strip()
                    context = (row.get('CONTEXT') or '').strip()
                    region = (row.get('REGION') or '').strip()
                    chr_id = (row.get('CHR_ID') or '').strip()
                    chr_pos = (row.get('CHR_POS') or '').strip()
                    risk_allele_freq = (
                        row.get('RISK ALLELE FREQUENCY') or ''
                    ).strip()
                    pvalue_mlog_str = (row.get('PVALUE_MLOG') or '').strip()

                    # Parse numeric properties
                    try:
                        or_beta = float(or_beta_str)
                    except (ValueError, TypeError):
                        or_beta = None

                    try:
                        risk_freq = float(risk_allele_freq)
                    except (ValueError, TypeError):
                        risk_freq = None

                    try:
                        pvalue_mlog = float(pvalue_mlog_str)
                    except (ValueError, TypeError):
                        pvalue_mlog = None

                    # Create one association per trait URI
                    for i, trait_uri in enumerate(trait_uris):
                        efo_id = self._extract_efo_id(trait_uri)
                        if not efo_id:
                            continue

                        trait_name = (
                            trait_names[i]
                            if i < len(trait_names)
                            else trait_names[0] if trait_names else ''
                        )

                        # Register trait node
                        if efo_id not in self.traits:
                            self.traits[efo_id] = {
                                'name': trait_name,
                                'uri': trait_uri,
                            }

                        assoc_data = {
                            'snp_id': snp_id,
                            'trait_id': efo_id,
                            'p_value': pval,
                            'pvalue_mlog': pvalue_mlog,
                            'or_beta': or_beta,
                            'confidence_interval': ci_text,
                            'risk_allele': risk_allele,
                            'risk_allele_frequency': risk_freq,
                            'pubmed_id': pubmed_id,
                            'study_accession': study_accession,
                            'mapped_gene': mapped_gene,
                            'context': context,
                            'region': region,
                            'chromosome': chr_id,
                            'position': chr_pos,
                        }

                        # Deduplicate: keep the most significant p-value
                        dedup_key = (snp_id, efo_id)
                        if dedup_key in best:
                            if pval < best[dedup_key]['p_value']:
                                best[dedup_key] = assoc_data
                        else:
                            best[dedup_key] = assoc_data

        self.associations = list(best.values())
        logger.info(
            f"GWAS Catalog: Loaded {len(self.associations)} significant "
            f"associations ({len(self.traits)} unique traits). "
            f"Skipped: {skipped_pval} non-significant, "
            f"{skipped_snp} missing/invalid SNP, "
            f"{skipped_trait} missing trait."
        )

    def get_nodes(self):
        """
        Generate Disease/Trait nodes from EFO ontology terms.
        Yields: (id, label, properties)
        """
        logger.info("GWAS Catalog: Generating trait nodes...")
        count = 0

        for efo_id, info in self.traits.items():
            props = {
                'name': self._sanitize(info['name']),
                'uri': self._sanitize(info['uri']),
                'source': 'GWAS Catalog',
            }

            yield (efo_id, "Disease", props)
            count += 1

        logger.info(f"GWAS Catalog: Generated {count} Disease/Trait nodes")

    def get_edges(self):
        """
        Generate GWASAssociation edges (variant to disease association).
        Yields: (id, source, target, label, properties)

        Source: dbsnp:<rs_id> (SNP/variant)
        Target: EFO:<id> (disease/trait)
        """
        logger.info("GWAS Catalog: Generating association edges...")
        count = 0

        for assoc in self.associations:
            source_id = f"dbsnp:{assoc['snp_id']}"
            target_id = assoc['trait_id']

            props = {
                'p_value': assoc['p_value'],
                'pvalue_mlog': (
                    assoc['pvalue_mlog']
                    if assoc['pvalue_mlog'] is not None
                    else 0.0
                ),
                'or_beta': (
                    assoc['or_beta']
                    if assoc['or_beta'] is not None
                    else 0.0
                ),
                'confidence_interval': self._sanitize(
                    assoc['confidence_interval']
                ),
                'risk_allele': self._sanitize(assoc['risk_allele']),
                'risk_allele_frequency': (
                    assoc['risk_allele_frequency']
                    if assoc['risk_allele_frequency'] is not None
                    else 0.0
                ),
                'pubmed_id': self._sanitize(assoc['pubmed_id']),
                'study_accession': self._sanitize(assoc['study_accession']),
                'mapped_gene': self._sanitize(
                    assoc['mapped_gene'][:200]
                    if assoc['mapped_gene']
                    else ''
                ),
                'context': self._sanitize(assoc['context']),
                'region': self._sanitize(assoc['region']),
                'chromosome': self._sanitize(assoc['chromosome']),
                'position': self._sanitize(assoc['position']),
                'source': 'GWAS Catalog',
            }

            yield (
                None,
                source_id,
                target_id,
                "variant to disease association",
                props,
            )
            count += 1

        logger.info(
            f"GWAS Catalog: Generated {count} "
            f"variant-disease association edges"
        )
