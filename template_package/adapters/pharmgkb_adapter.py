"""
PharmGKB Adapter for BioCypher.

Loads pharmacogenomics data linking genes, drugs, and diseases.
Generates:
- Drug nodes from drugs.tsv (3.7K drugs)
- Gene-Drug association edges from relationships.tsv
- Gene-Disease association edges from relationships.tsv

PharmGKB is a comprehensive pharmacogenomics knowledge resource
that curates knowledge about the impact of genetic variation on
drug response for clinicians and researchers.
"""

import csv
import io
import zipfile
from pathlib import Path
from biocypher._logger import logger

# PharmGKB can have very large fields (e.g. cross-references, InChI strings)
csv.field_size_limit(10 * 1024 * 1024)


class PharmGKBAdapter:
    def __init__(self, data_dir="template_package/data/pharmgkb"):
        self.data_dir = Path(data_dir)
        self.drugs = []
        self.gene_drug_edges = []
        self.gene_disease_edges = []
        self._load_data()

    # ------------------------------------------------------------------
    # Sanitisation
    # ------------------------------------------------------------------

    def _sanitize(self, text):
        """Clean a text value for safe downstream use in BioCypher."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
        return text.strip()

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _read_tsv_from_zip(self, zip_filename, tsv_filename):
        """Open a TSV file inside a zip archive and return a csv.DictReader."""
        zip_path = self.data_dir / zip_filename
        if not zip_path.exists():
            logger.warning(f"PharmGKB: {zip_filename} not found at {zip_path}")
            return None, None
        zf = zipfile.ZipFile(zip_path, "r")
        raw = zf.open(tsv_filename)
        text_wrapper = io.TextIOWrapper(raw, encoding="utf-8")
        reader = csv.DictReader(text_wrapper, delimiter="\t")
        return reader, zf

    def _load_data(self):
        """Load drugs, genes (for reference), and relationships."""
        self._load_drugs()
        self._load_relationships()

    def _load_drugs(self):
        """Load drug nodes from drugs.zip / drugs.tsv."""
        reader, zf = self._read_tsv_from_zip("drugs.zip", "drugs.tsv")
        if reader is None:
            return

        logger.info("PharmGKB: Loading drugs...")
        count = 0

        try:
            for row in reader:
                drug_id = (row.get("PharmGKB Accession Id") or "").strip()
                if not drug_id:
                    continue

                name = (row.get("Name") or "").strip()
                generic_names = (row.get("Generic Names") or "").strip()
                trade_names = (row.get("Trade Names") or "").strip()
                brand_mixtures = (row.get("Brand Mixtures") or "").strip()
                drug_type = (row.get("Type") or "").strip()
                cross_refs = (row.get("Cross-references") or "").strip()
                smiles = (row.get("SMILES") or "").strip()
                inchi = (row.get("InChI") or "").strip()
                dosing_guideline = (row.get("Dosing Guideline") or "").strip()
                external_vocab = (row.get("External Vocabulary") or "").strip()
                clinical_ann_count = (row.get("Clinical Annotation Count") or "").strip()
                variant_ann_count = (row.get("Variant Annotation Count") or "").strip()
                pathway_count = (row.get("Pathway Count") or "").strip()
                atc_ids = (row.get("ATC Identifiers") or "").strip()
                rxnorm_ids = (row.get("RxNorm Identifiers") or "").strip()
                pubchem_ids = (row.get("PubChem Compound Identifiers") or "").strip()
                top_clin_level = (row.get("Top Clinical Annotation Level") or "").strip()
                top_fda_level = (row.get("Top FDA Label Testing Level") or "").strip()
                top_cpic_level = (row.get("Top CPIC Pairs Level") or "").strip()
                dosing_sources = (row.get("Dosing Guideline Sources") or "").strip()

                self.drugs.append({
                    "drug_id": drug_id,
                    "name": name,
                    "generic_names": generic_names,
                    "trade_names": trade_names,
                    "brand_mixtures": brand_mixtures,
                    "drug_type": drug_type,
                    "cross_references": cross_refs[:500],
                    "smiles": smiles,
                    "inchi": inchi[:500],
                    "dosing_guideline": dosing_guideline,
                    "external_vocabulary": external_vocab[:300],
                    "clinical_annotation_count": clinical_ann_count,
                    "variant_annotation_count": variant_ann_count,
                    "pathway_count": pathway_count,
                    "atc_identifiers": atc_ids,
                    "rxnorm_identifiers": rxnorm_ids,
                    "pubchem_compound_ids": pubchem_ids,
                    "top_clinical_annotation_level": top_clin_level,
                    "top_fda_label_testing_level": top_fda_level,
                    "top_cpic_pairs_level": top_cpic_level,
                    "dosing_guideline_sources": dosing_sources,
                })
                count += 1
        finally:
            zf.close()

        logger.info(f"PharmGKB: Loaded {count} drugs")

    def _load_relationships(self):
        """
        Load gene-drug and gene-disease relationship edges from
        relationships.zip / relationships.tsv.

        The file contains relationships between many entity type pairs
        (Gene, Chemical, Disease, Variant, Haplotype). We extract:
        - Gene <-> Chemical  => gene_drug_edges
        - Gene <-> Disease   => gene_disease_edges
        """
        reader, zf = self._read_tsv_from_zip("relationships.zip", "relationships.tsv")
        if reader is None:
            return

        logger.info("PharmGKB: Loading relationships...")
        gene_drug_count = 0
        gene_disease_count = 0
        seen_gene_drug = set()
        seen_gene_disease = set()

        try:
            for row in reader:
                e1_id = (row.get("Entity1_id") or "").strip()
                e1_name = (row.get("Entity1_name") or "").strip()
                e1_type = (row.get("Entity1_type") or "").strip()
                e2_id = (row.get("Entity2_id") or "").strip()
                e2_name = (row.get("Entity2_name") or "").strip()
                e2_type = (row.get("Entity2_type") or "").strip()
                evidence = (row.get("Evidence") or "").strip()
                association = (row.get("Association") or "").strip()
                pk = (row.get("PK") or "").strip()
                pd = (row.get("PD") or "").strip()
                pmids = (row.get("PMIDs") or "").strip()

                if not e1_id or not e2_id:
                    continue

                # Normalise so gene is always the source entity
                # Gene <-> Chemical (Drug)
                if e1_type == "Gene" and e2_type == "Chemical":
                    gene_id, gene_name = e1_id, e1_name
                    drug_id, drug_name = e2_id, e2_name
                elif e1_type == "Chemical" and e2_type == "Gene":
                    gene_id, gene_name = e2_id, e2_name
                    drug_id, drug_name = e1_id, e1_name
                # Gene <-> Disease
                elif e1_type == "Gene" and e2_type == "Disease":
                    gene_id, gene_name = e1_id, e1_name
                    disease_id, disease_name = e2_id, e2_name
                elif e1_type == "Disease" and e2_type == "Gene":
                    gene_id, gene_name = e2_id, e2_name
                    disease_id, disease_name = e1_id, e1_name
                else:
                    # Skip other entity type pairs (Variant, Haplotype, etc.)
                    continue

                # ---- Gene-Drug edges ----
                if (e1_type == "Gene" and e2_type == "Chemical") or \
                   (e1_type == "Chemical" and e2_type == "Gene"):
                    dedup_key = (gene_id, drug_id)
                    if dedup_key in seen_gene_drug:
                        continue
                    seen_gene_drug.add(dedup_key)

                    self.gene_drug_edges.append({
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "drug_id": drug_id,
                        "drug_name": drug_name,
                        "evidence": evidence,
                        "association": association,
                        "pk": pk,
                        "pd": pd,
                        "pmids": pmids,
                    })
                    gene_drug_count += 1

                # ---- Gene-Disease edges ----
                elif (e1_type == "Gene" and e2_type == "Disease") or \
                     (e1_type == "Disease" and e2_type == "Gene"):
                    dedup_key = (gene_id, disease_id)
                    if dedup_key in seen_gene_disease:
                        continue
                    seen_gene_disease.add(dedup_key)

                    self.gene_disease_edges.append({
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "disease_id": disease_id,
                        "disease_name": disease_name,
                        "evidence": evidence,
                        "association": association,
                        "pk": pk,
                        "pd": pd,
                        "pmids": pmids,
                    })
                    gene_disease_count += 1
        finally:
            zf.close()

        logger.info(
            f"PharmGKB: Loaded {gene_drug_count} gene-drug and "
            f"{gene_disease_count} gene-disease relationships"
        )

    # ------------------------------------------------------------------
    # Node generation
    # ------------------------------------------------------------------

    def get_nodes(self):
        """
        Generate Drug nodes.

        Yields: (id, label, properties)
            - id:    "pharmgkb:<PharmGKB Accession Id>"
            - label: "drug"
            - properties: dict of drug attributes
        """
        logger.info("PharmGKB: Generating Drug nodes...")
        count = 0

        for drug in self.drugs:
            props = {
                "name": self._sanitize(drug["name"]),
                "generic_names": self._sanitize(drug["generic_names"]),
                "trade_names": self._sanitize(drug["trade_names"]),
                "brand_mixtures": self._sanitize(drug["brand_mixtures"]),
                "drug_type": self._sanitize(drug["drug_type"]),
                "cross_references": self._sanitize(drug["cross_references"]),
                "smiles": self._sanitize(drug["smiles"]),
                "inchi": self._sanitize(drug["inchi"]),
                "dosing_guideline": self._sanitize(drug["dosing_guideline"]),
                "external_vocabulary": self._sanitize(drug["external_vocabulary"]),
                "clinical_annotation_count": drug["clinical_annotation_count"],
                "variant_annotation_count": drug["variant_annotation_count"],
                "pathway_count": drug["pathway_count"],
                "atc_identifiers": self._sanitize(drug["atc_identifiers"]),
                "rxnorm_identifiers": self._sanitize(drug["rxnorm_identifiers"]),
                "pubchem_compound_ids": self._sanitize(drug["pubchem_compound_ids"]),
                "top_clinical_annotation_level": drug["top_clinical_annotation_level"],
                "top_fda_label_testing_level": drug["top_fda_label_testing_level"],
                "top_cpic_pairs_level": drug["top_cpic_pairs_level"],
                "dosing_guideline_sources": self._sanitize(drug["dosing_guideline_sources"]),
                "source": "PharmGKB",
            }

            yield (f"pharmgkb:{drug['drug_id']}", "drug", props)
            count += 1

        logger.info(f"PharmGKB: Generated {count} Drug nodes")

    # ------------------------------------------------------------------
    # Edge generation
    # ------------------------------------------------------------------

    def get_edges(self):
        """
        Generate gene-drug and gene-disease association edges.

        Yields: (id, source, target, label, properties)
            - id:     None (auto-generated)
            - source: "pharmgkb:<gene_id>"
            - target: "pharmgkb:<drug_id>" or "pharmgkb:<disease_id>"
            - label:  "gene to drug association" or "gene to disease association"
            - properties: evidence, association, pk_pd, pubmed_ids, names
        """
        logger.info("PharmGKB: Generating edges...")
        count = 0

        # Gene-Drug associations
        for edge in self.gene_drug_edges:
            pk_pd_parts = []
            if edge["pk"]:
                pk_pd_parts.append(f"PK:{edge['pk']}")
            if edge["pd"]:
                pk_pd_parts.append(f"PD:{edge['pd']}")
            pk_pd = "; ".join(pk_pd_parts)

            props = {
                "gene_name": self._sanitize(edge["gene_name"]),
                "drug_name": self._sanitize(edge["drug_name"]),
                "evidence": self._sanitize(edge["evidence"]),
                "association": self._sanitize(edge["association"]),
                "pk_pd": self._sanitize(pk_pd),
                "pubmed_ids": self._sanitize(edge["pmids"]),
                "source": "PharmGKB",
            }

            yield (
                None,
                f"pharmgkb:{edge['gene_id']}",
                f"pharmgkb:{edge['drug_id']}",
                "gene to drug association",
                props,
            )
            count += 1

        gene_drug_count = count

        # Gene-Disease associations
        for edge in self.gene_disease_edges:
            pk_pd_parts = []
            if edge["pk"]:
                pk_pd_parts.append(f"PK:{edge['pk']}")
            if edge["pd"]:
                pk_pd_parts.append(f"PD:{edge['pd']}")
            pk_pd = "; ".join(pk_pd_parts)

            props = {
                "gene_name": self._sanitize(edge["gene_name"]),
                "disease_name": self._sanitize(edge["disease_name"]),
                "evidence": self._sanitize(edge["evidence"]),
                "association": self._sanitize(edge["association"]),
                "pk_pd": self._sanitize(pk_pd),
                "pubmed_ids": self._sanitize(edge["pmids"]),
                "source": "PharmGKB",
            }

            yield (
                None,
                f"pharmgkb:{edge['gene_id']}",
                f"pharmgkb:{edge['disease_id']}",
                "gene to disease association",
                props,
            )
            count += 1

        gene_disease_count = count - gene_drug_count

        logger.info(
            f"PharmGKB: Generated {gene_drug_count} gene-drug and "
            f"{gene_disease_count} gene-disease edges ({count} total)"
        )
