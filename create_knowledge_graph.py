"""
Main script to create the unified knowledge graph.
Loads all adapters and writes nodes/edges through BioCypher.
"""

import biocypher
from biocypher._logger import logger

# ============================================================
# 1. Instantiate BioCypher Driver
# ============================================================
driver = biocypher.BioCypher()

# ============================================================
# 2. Define adapters to load (in order of ontology dependency)
# ============================================================
adapters = []

# --- LIANA (Ligand-Receptor Interactions) ---
try:
    from template_package.adapters.liana_adapter import LianaAdapter
    adapters.append(("LIANA", LianaAdapter()))
    logger.info("Loaded LIANA adapter")
except Exception as e:
    logger.warning(f"Could not load LIANA adapter: {e}")

# --- ChEBI (Chemical Entities) ---
# Foundational: provides chemical vocabulary for other databases
try:
    from template_package.adapters.chebi_adapter import ChEBIAdapter
    adapters.append(("ChEBI", ChEBIAdapter()))
    logger.info("Loaded ChEBI adapter")
except Exception as e:
    logger.warning(f"Could not load ChEBI adapter: {e}")

# --- Rhea (Biochemical Reactions) ---
# Depends on ChEBI for substrate/product linking
try:
    from template_package.adapters.rhea_adapter import RheaAdapter
    adapters.append(("Rhea", RheaAdapter()))
    logger.info("Loaded Rhea adapter")
except Exception as e:
    logger.warning(f"Could not load Rhea adapter: {e}")

# --- LIPID MAPS (Lipid Structures) ---
# Cross-references to ChEBI
try:
    from template_package.adapters.lipidmaps_adapter import LIPIDMAPSAdapter
    adapters.append(("LIPIDMAPS", LIPIDMAPSAdapter()))
    logger.info("Loaded LIPIDMAPS adapter")
except Exception as e:
    logger.warning(f"Could not load LIPIDMAPS adapter: {e}")

# --- KEGG (Pathways, Reactions, Compounds) ---
# Cross-references to ChEBI
try:
    from template_package.adapters.kegg_adapter import KEGGAdapter
    adapters.append(("KEGG", KEGGAdapter()))
    logger.info("Loaded KEGG adapter")
except Exception as e:
    logger.warning(f"Could not load KEGG adapter: {e}")

# --- ComplexPortal (Protein Complexes) ---
try:
    from template_package.adapters.complexportal_adapter import ComplexPortalAdapter
    adapters.append(("ComplexPortal", ComplexPortalAdapter()))
    logger.info("Loaded ComplexPortal adapter")
except Exception as e:
    logger.warning(f"Could not load ComplexPortal adapter: {e}")

# --- InterPro (Protein Families and Domains) ---
try:
    from template_package.adapters.interpro_adapter import InterProAdapter
    adapters.append(("InterPro", InterProAdapter()))
    logger.info("Loaded InterPro adapter")
except Exception as e:
    logger.warning(f"Could not load InterPro adapter: {e}")

# --- STRING (Protein-Protein Interactions) ---
try:
    from template_package.adapters.string_adapter import STRINGAdapter
    adapters.append(("STRING", STRINGAdapter()))
    logger.info("Loaded STRING adapter")
except Exception as e:
    logger.warning(f"Could not load STRING adapter: {e}")

# --- Reactome (Biological Pathways) ---
try:
    from template_package.adapters.reactome_adapter import ReactomeAdapter
    adapters.append(("Reactome", ReactomeAdapter()))
    logger.info("Loaded Reactome adapter")
except Exception as e:
    logger.warning(f"Could not load Reactome adapter: {e}")

# --- Human Protein Atlas (Expression & Localization) ---
try:
    from template_package.adapters.hpa_adapter import HPAAdapter
    adapters.append(("HPA", HPAAdapter()))
    logger.info("Loaded HPA adapter")
except Exception as e:
    logger.warning(f"Could not load HPA adapter: {e}")

# --- IntAct (Molecular Interactions) ---
try:
    from template_package.adapters.intact_adapter import IntActAdapter
    adapters.append(("IntAct", IntActAdapter()))
    logger.info("Loaded IntAct adapter")
except Exception as e:
    logger.warning(f"Could not load IntAct adapter: {e}")

# --- Compartments (Jensen Lab - Subcellular Localization) ---
try:
    from template_package.adapters.compartments_adapter import CompartmentsAdapter
    adapters.append(("Compartments", CompartmentsAdapter()))
    logger.info("Loaded Compartments adapter")
except Exception as e:
    logger.warning(f"Could not load Compartments adapter: {e}")

# --- GPCRdb (G Protein-Coupled Receptors) ---
try:
    from template_package.adapters.gpcrdb_adapter import GPCRdbAdapter
    adapters.append(("GPCRdb", GPCRdbAdapter()))
    logger.info("Loaded GPCRdb adapter")
except Exception as e:
    logger.warning(f"Could not load GPCRdb adapter: {e}")

# --- ELM (Eukaryotic Linear Motifs) ---
try:
    from template_package.adapters.elm_adapter import ELMAdapter
    adapters.append(("ELM", ELMAdapter()))
    logger.info("Loaded ELM adapter")
except Exception as e:
    logger.warning(f"Could not load ELM adapter: {e}")

# --- Cell Ontology ---
try:
    from template_package.adapters.cell_ontology_adapter import CellOntologyAdapter
    adapters.append(("CellOntology", CellOntologyAdapter()))
    logger.info("Loaded CellOntology adapter")
except Exception as e:
    logger.warning(f"Could not load CellOntology adapter: {e}")

# --- Allen Brain Atlas (Brain Anatomy) ---
try:
    from template_package.adapters.allen_brain_adapter import AllenBrainAdapter
    adapters.append(("AllenBrain", AllenBrainAdapter()))
    logger.info("Loaded AllenBrain adapter")
except Exception as e:
    logger.warning(f"Could not load AllenBrain adapter: {e}")

# --- NeuroElectro (Neuronal Electrophysiology) ---
try:
    from template_package.adapters.neuroelectro_adapter import NeuroElectroAdapter
    adapters.append(("NeuroElectro", NeuroElectroAdapter()))
    logger.info("Loaded NeuroElectro adapter")
except Exception as e:
    logger.warning(f"Could not load NeuroElectro adapter: {e}")

# --- SynGO (Synaptic Gene Ontologies) ---
try:
    from template_package.adapters.syngo_adapter import SynGOAdapter
    adapters.append(("SynGO", SynGOAdapter()))
    logger.info("Loaded SynGO adapter")
except Exception as e:
    logger.warning(f"Could not load SynGO adapter: {e}")

# --- ENCODE SCREEN (Candidate cis-Regulatory Elements) ---
try:
    from template_package.adapters.encode_screen_adapter import ENCODESCREENAdapter
    adapters.append(("ENCODE_SCREEN", ENCODESCREENAdapter()))
    logger.info("Loaded ENCODE SCREEN adapter")
except Exception as e:
    logger.warning(f"Could not load ENCODE SCREEN adapter: {e}")

# --- BioLiP (Ligand-Protein Interactions) ---
try:
    from template_package.adapters.biolip_adapter import BioLiPAdapter
    adapters.append(("BioLiP", BioLiPAdapter()))
    logger.info("Loaded BioLiP adapter")
except Exception as e:
    logger.warning(f"Could not load BioLiP adapter: {e}")

# --- TCDB (Transporter Classification Database) ---
try:
    from template_package.adapters.tcdb_adapter import TCDBAdapter
    adapters.append(("TCDB", TCDBAdapter()))
    logger.info("Loaded TCDB adapter")
except Exception as e:
    logger.warning(f"Could not load TCDB adapter: {e}")

# --- DGIdb (Drug Gene Interaction Database) ---
try:
    from template_package.adapters.dgidb_adapter import DGIdbAdapter
    adapters.append(("DGIdb", DGIdbAdapter()))
    logger.info("Loaded DGIdb adapter")
except Exception as e:
    logger.warning(f"Could not load DGIdb adapter: {e}")

# --- GlyGen (Glycan Structures) ---
try:
    from template_package.adapters.glygen_adapter import GlyGenAdapter
    adapters.append(("GlyGen", GlyGenAdapter()))
    logger.info("Loaded GlyGen adapter")
except Exception as e:
    logger.warning(f"Could not load GlyGen adapter: {e}")

# --- MatrisomeDB (Extracellular Matrix) ---
try:
    from template_package.adapters.matrisome_adapter import MatrisomeAdapter
    adapters.append(("MatrisomeDB", MatrisomeAdapter()))
    logger.info("Loaded MatrisomeDB adapter")
except Exception as e:
    logger.warning(f"Could not load MatrisomeDB adapter: {e}")

# --- ExoCarta (Exosome Database) ---
try:
    from template_package.adapters.exocarta_adapter import ExoCartaAdapter
    adapters.append(("ExoCarta", ExoCartaAdapter()))
    logger.info("Loaded ExoCarta adapter")
except Exception as e:
    logger.warning(f"Could not load ExoCarta adapter: {e}")

# --- CPDB (ConsensusPathDB) ---
try:
    from template_package.adapters.cpdb_adapter import CPDBAdapter
    adapters.append(("CPDB", CPDBAdapter()))
    logger.info("Loaded CPDB adapter")
except Exception as e:
    logger.warning(f"Could not load CPDB adapter: {e}")

# --- MobiDB (Intrinsic Disorder) ---
try:
    from template_package.adapters.mobidb_adapter import MobiDBAdapter
    adapters.append(("MobiDB", MobiDBAdapter()))
    logger.info("Loaded MobiDB adapter")
except Exception as e:
    logger.warning(f"Could not load MobiDB adapter: {e}")

# --- iPTMnet (Post-Translational Modifications) ---
try:
    from template_package.adapters.iptmnet_adapter import iPTMnetAdapter
    adapters.append(("iPTMnet", iPTMnetAdapter()))
    logger.info("Loaded iPTMnet adapter")
except Exception as e:
    logger.warning(f"Could not load iPTMnet adapter: {e}")

# --- Dfam (Transposable Element Families) ---
try:
    from template_package.adapters.dfam_adapter import DfamAdapter
    adapters.append(("Dfam", DfamAdapter()))
    logger.info("Loaded Dfam adapter")
except Exception as e:
    logger.warning(f"Could not load Dfam adapter: {e}")

# --- REPAIRtoire (DNA Repair Pathways) ---
try:
    from template_package.adapters.repairtoire_adapter import RepairToireAdapter
    adapters.append(("REPAIRtoire", RepairToireAdapter()))
    logger.info("Loaded REPAIRtoire adapter")
except Exception as e:
    logger.warning(f"Could not load REPAIRtoire adapter: {e}")

# --- ComPPI (Compartmentalized PPI) ---
try:
    from template_package.adapters.comppi_adapter import ComPPIAdapter
    adapters.append(("ComPPI", ComPPIAdapter()))
    logger.info("Loaded ComPPI adapter")
except Exception as e:
    logger.warning(f"Could not load ComPPI adapter: {e}")

# --- FANTOM5 (Enhancers) ---
try:
    from template_package.adapters.fantom5_adapter import FANTOM5Adapter
    adapters.append(("FANTOM5", FANTOM5Adapter()))
    logger.info("Loaded FANTOM5 adapter")
except Exception as e:
    logger.warning(f"Could not load FANTOM5 adapter: {e}")

# --- Membranome (Transmembrane Proteins) ---
try:
    from template_package.adapters.membranome_adapter import MembranomeAdapter
    adapters.append(("Membranome", MembranomeAdapter()))
    logger.info("Loaded Membranome adapter")
except Exception as e:
    logger.warning(f"Could not load Membranome adapter: {e}")

# --- SLC Tables (Solute Carriers) ---
try:
    from template_package.adapters.slc_adapter import SLCAdapter
    adapters.append(("SLC", SLCAdapter()))
    logger.info("Loaded SLC adapter")
except Exception as e:
    logger.warning(f"Could not load SLC adapter: {e}")

# --- HistoneDB (Histone Variants) ---
try:
    from template_package.adapters.histonedb_adapter import HistoneDBAdapter
    adapters.append(("HistoneDB", HistoneDBAdapter()))
    logger.info("Loaded HistoneDB adapter")
except Exception as e:
    logger.warning(f"Could not load HistoneDB adapter: {e}")

# --- MetalPDB (Metal Binding Sites) ---
try:
    from template_package.adapters.metalpdb_adapter import MetalPDBAdapter
    adapters.append(("MetalPDB", MetalPDBAdapter()))
    logger.info("Loaded MetalPDB adapter")
except Exception as e:
    logger.warning(f"Could not load MetalPDB adapter: {e}")

# --- ChEA3 (TF Target Interactions) ---
try:
    from template_package.adapters.chea3_adapter import ChEA3Adapter
    adapters.append(("ChEA3", ChEA3Adapter()))
    logger.info("Loaded ChEA3 adapter")
except Exception as e:
    logger.warning(f"Could not load ChEA3 adapter: {e}")

# --- OpenProt (Alternative Proteins) ---
try:
    from template_package.adapters.openprot_adapter import OpenProtAdapter
    adapters.append(("OpenProt", OpenProtAdapter()))
    logger.info("Loaded OpenProt adapter")
except Exception as e:
    logger.warning(f"Could not load OpenProt adapter: {e}")

# --- PTMcode (PTM Crosstalk) ---
try:
    from template_package.adapters.ptmcode_adapter import PTMcodeAdapter
    adapters.append(("PTMcode", PTMcodeAdapter()))
    logger.info("Loaded PTMcode adapter")
except Exception as e:
    logger.warning(f"Could not load PTMcode adapter: {e}")

# --- ModelDB (Computational Neuroscience Models) ---
try:
    from template_package.adapters.modeldb_adapter import ModelDBAdapter
    adapters.append(("ModelDB", ModelDBAdapter()))
    logger.info("Loaded ModelDB adapter")
except Exception as e:
    logger.warning(f"Could not load ModelDB adapter: {e}")

# --- PeptideAtlas (Proteomics Detection) ---
try:
    from template_package.adapters.peptideatlas_adapter import PeptideAtlasAdapter
    adapters.append(("PeptideAtlas", PeptideAtlasAdapter()))
    logger.info("Loaded PeptideAtlas adapter")
except Exception as e:
    logger.warning(f"Could not load PeptideAtlas adapter: {e}")

# --- SABIO-RK (Enzyme Kinetics) ---
try:
    from template_package.adapters.sabiork_adapter import SABIORKAdapter
    adapters.append(("SABIORK", SABIORKAdapter()))
    logger.info("Loaded SABIO-RK adapter")
except Exception as e:
    logger.warning(f"Could not load SABIO-RK adapter: {e}")

# --- PDBe PISA (Protein Interfaces) ---
try:
    from template_package.adapters.pdbe_pisa_adapter import PDBePISAAdapter
    adapters.append(("PDBePISA", PDBePISAAdapter()))
    logger.info("Loaded PDBe PISA adapter")
except Exception as e:
    logger.warning(f"Could not load PDBe PISA adapter: {e}")

# --- ExPASy ENZYME (Enzyme Nomenclature) ---
try:
    from template_package.adapters.enzyme_adapter import EnzymeAdapter
    adapters.append(("ENZYME", EnzymeAdapter()))
    logger.info("Loaded ExPASy ENZYME adapter")
except Exception as e:
    logger.warning(f"Could not load ENZYME adapter: {e}")

# --- PhaSepDB (Phase Separation) ---
try:
    from template_package.adapters.phasep_adapter import PhasePAdapter
    adapters.append(("PhaseP", PhasePAdapter()))
    logger.info("Loaded PhaSepDB adapter")
except Exception as e:
    logger.warning(f"Could not load PhaseP adapter: {e}")

# --- RNAgranuleDB (RNA Granules) ---
try:
    from template_package.adapters.rnagranuledb_adapter import RNAgranuleDBAdapter
    adapters.append(("RNAgranuleDB", RNAgranuleDBAdapter()))
    logger.info("Loaded RNAgranuleDB adapter")
except Exception as e:
    logger.warning(f"Could not load RNAgranuleDB adapter: {e}")

# --- Degronopedia (Degron Motifs) ---
try:
    from template_package.adapters.degronopedia_adapter import DegronopediaAdapter
    adapters.append(("Degronopedia", DegronopediaAdapter()))
    logger.info("Loaded Degronopedia adapter")
except Exception as e:
    logger.warning(f"Could not load Degronopedia adapter: {e}")

# --- FerrDb (Ferroptosis Genes) ---
try:
    from template_package.adapters.ferrdb_adapter import FerrDbAdapter
    adapters.append(("FerrDb", FerrDbAdapter()))
    logger.info("Loaded FerrDb adapter")
except Exception as e:
    logger.warning(f"Could not load FerrDb adapter: {e}")

# --- PsychENCODE (Brain Regulatory Elements) ---
try:
    from template_package.adapters.psychencode_adapter import PsychENCODEAdapter
    adapters.append(("PsychENCODE", PsychENCODEAdapter()))
    logger.info("Loaded PsychENCODE adapter")
except Exception as e:
    logger.warning(f"Could not load PsychENCODE adapter: {e}")

# --- ArrestinDB (GPCR-Arrestin Coupling) ---
try:
    from template_package.adapters.arrestindb_adapter import ArrestinDBAdapter
    adapters.append(("ArrestinDB", ArrestinDBAdapter()))
    logger.info("Loaded ArrestinDB adapter")
except Exception as e:
    logger.warning(f"Could not load ArrestinDB adapter: {e}")

# --- GtRNAdb (tRNA Genes) ---
try:
    from template_package.adapters.gtrnadb_adapter import GtRNAdbAdapter
    adapters.append(("GtRNAdb", GtRNAdbAdapter()))
    logger.info("Loaded GtRNAdb adapter")
except Exception as e:
    logger.warning(f"Could not load GtRNAdb adapter: {e}")

# --- Allen Brain Connectivity Atlas ---
try:
    from template_package.adapters.allen_connectivity_adapter import AllenConnectivityAdapter
    adapters.append(("AllenConnectivity", AllenConnectivityAdapter()))
    logger.info("Loaded Allen Connectivity adapter")
except Exception as e:
    logger.warning(f"Could not load Allen Connectivity adapter: {e}")

# --- Translocatome (Protein Translocation) ---
try:
    from template_package.adapters.translocatome_adapter import TranslocatomeAdapter
    adapters.append(("Translocatome", TranslocatomeAdapter()))
    logger.info("Loaded Translocatome adapter")
except Exception as e:
    logger.warning(f"Could not load Translocatome adapter: {e}")


# --- 3D Genome Browser (TADs) ---
try:
    from template_package.adapters.threed_genome_adapter import ThreeDGenomeAdapter
    adapters.append(("3DGenome", ThreeDGenomeAdapter()))
    logger.info("Loaded 3D Genome Browser adapter")
except Exception as e:
    logger.warning(f"Could not load 3D Genome Browser adapter: {e}")

# --- Adhesome (Cell Adhesion) ---
try:
    from template_package.adapters.adhesome_adapter import AdhesomeAdapter
    adapters.append(("Adhesome", AdhesomeAdapter()))
    logger.info("Loaded Adhesome adapter")
except Exception as e:
    logger.warning(f"Could not load Adhesome adapter: {e}")

# --- BRENDA (Enzyme Database) ---
try:
    from template_package.adapters.brenda_adapter import BRENDAAdapter
    adapters.append(("BRENDA", BRENDAAdapter()))
    logger.info("Loaded BRENDA adapter")
except Exception as e:
    logger.warning(f"Could not load BRENDA adapter: {e}")

# --- ChapNet (Chaperone Interaction Network) ---
try:
    from template_package.adapters.chapnet_adapter import ChapNetAdapter
    adapters.append(("ChapNet", ChapNetAdapter()))
    logger.info("Loaded ChapNet adapter")
except Exception as e:
    logger.warning(f"Could not load ChapNet adapter: {e}")

# --- ENCORI / starBase (miRNA-Target Interactions) ---
try:
    from template_package.adapters.encori_adapter import ENCORIAdapter
    adapters.append(("ENCORI", ENCORIAdapter()))
    logger.info("Loaded ENCORI adapter")
except Exception as e:
    logger.warning(f"Could not load ENCORI adapter: {e}")

# --- EpiMap (Chromatin State Segmentation) ---
try:
    from template_package.adapters.epimap_adapter import EpiMapAdapter
    adapters.append(("EpiMap", EpiMapAdapter()))
    logger.info("Loaded EpiMap adapter")
except Exception as e:
    logger.warning(f"Could not load EpiMap adapter: {e}")

# --- EVpedia (Extracellular Vesicle) ---
try:
    from template_package.adapters.evpedia_adapter import EVpediaAdapter
    adapters.append(("EVpedia", EVpediaAdapter()))
    logger.info("Loaded EVpedia adapter")
except Exception as e:
    logger.warning(f"Could not load EVpedia adapter: {e}")

# --- LION (Lipid Ontology) ---
try:
    from template_package.adapters.lipid_ontology_adapter import LipidOntologyAdapter
    adapters.append(("LION", LipidOntologyAdapter()))
    logger.info("Loaded Lipid Ontology (LION) adapter")
except Exception as e:
    logger.warning(f"Could not load Lipid Ontology adapter: {e}")

# --- RaftProt (Lipid Raft Proteomics) ---
try:
    from template_package.adapters.raftprot_adapter import RaftProtAdapter
    adapters.append(("RaftProt", RaftProtAdapter()))
    logger.info("Loaded RaftProt adapter")
except Exception as e:
    logger.warning(f"Could not load RaftProt adapter: {e}")

# --- SILVA (rRNA Taxonomy) ---
try:
    from template_package.adapters.silva_adapter import SILVAAdapter
    adapters.append(("SILVA", SILVAAdapter()))
    logger.info("Loaded SILVA adapter")
except Exception as e:
    logger.warning(f"Could not load SILVA adapter: {e}")

# --- TarBase / miRDB (miRNA Target Predictions) ---
try:
    from template_package.adapters.tarbase_adapter import TarBaseAdapter
    adapters.append(("TarBase", TarBaseAdapter()))
    logger.info("Loaded TarBase/miRDB adapter")
except Exception as e:
    logger.warning(f"Could not load TarBase/miRDB adapter: {e}")

# --- UniCarbKB / GlyGen (Glycan Structures) ---
try:
    from template_package.adapters.unicarbkb_adapter import UniCarbKBAdapter
    adapters.append(("UniCarbKB", UniCarbKBAdapter()))
    logger.info("Loaded UniCarbKB/GlyGen adapter")
except Exception as e:
    logger.warning(f"Could not load UniCarbKB/GlyGen adapter: {e}")

# --- GENCODE (Gene Annotations) ---
try:
    from template_package.adapters.gencode_adapter import GENCODEAdapter
    adapters.append(("GENCODE", GENCODEAdapter()))
    logger.info("Loaded GENCODE adapter")
except Exception as e:
    logger.warning(f"Could not load GENCODE adapter: {e}")

# --- JASPAR (TF Binding Profiles) ---
try:
    from template_package.adapters.jaspar_adapter import JASPARAdapter
    adapters.append(("JASPAR", JASPARAdapter()))
    logger.info("Loaded JASPAR adapter")
except Exception as e:
    logger.warning(f"Could not load JASPAR adapter: {e}")

# --- HuRI (Human Reference Interactome) ---
try:
    from template_package.adapters.huri_adapter import HuRIAdapter
    adapters.append(("HuRI", HuRIAdapter()))
    logger.info("Loaded HuRI adapter")
except Exception as e:
    logger.warning(f"Could not load HuRI adapter: {e}")

# --- APPRIS (Principal Isoforms) ---
try:
    from template_package.adapters.appris_adapter import APPRISAdapter
    adapters.append(("APPRIS", APPRISAdapter()))
    logger.info("Loaded APPRIS adapter")
except Exception as e:
    logger.warning(f"Could not load APPRIS adapter: {e}")

# --- RNAcentral (Non-coding RNAs) ---
try:
    from template_package.adapters.rnacentral_adapter import RNAcentralAdapter
    adapters.append(("RNAcentral", RNAcentralAdapter()))
    logger.info("Loaded RNAcentral adapter")
except Exception as e:
    logger.warning(f"Could not load RNAcentral adapter: {e}")

# --- MODOMICS (RNA Modifications) ---
try:
    from template_package.adapters.modomics_adapter import MODOMICSAdapter
    adapters.append(("MODOMICS", MODOMICSAdapter()))
    logger.info("Loaded MODOMICS adapter")
except Exception as e:
    logger.warning(f"Could not load MODOMICS adapter: {e}")

# --- OPM (Membrane Protein Orientations) ---
try:
    from template_package.adapters.opm_adapter import OPMAdapter
    adapters.append(("OPM", OPMAdapter()))
    logger.info("Loaded OPM adapter")
except Exception as e:
    logger.warning(f"Could not load OPM adapter: {e}")

# --- PolyASite (Polyadenylation Sites) ---
try:
    from template_package.adapters.polyasite_adapter import PolyASiteAdapter
    adapters.append(("PolyASite", PolyASiteAdapter()))
    logger.info("Loaded PolyASite adapter")
except Exception as e:
    logger.warning(f"Could not load PolyASite adapter: {e}")

# --- DrLLPS (Liquid-Liquid Phase Separation) ---
try:
    from template_package.adapters.drllps_adapter import DrLLPSAdapter
    adapters.append(("DrLLPS", DrLLPSAdapter()))
    logger.info("Loaded DrLLPS adapter")
except Exception as e:
    logger.warning(f"Could not load DrLLPS adapter: {e}")

# --- ATtRACT (RNA-Binding Protein Motifs) ---
try:
    from template_package.adapters.attract_adapter import ATtRACTAdapter
    adapters.append(("ATtRACT", ATtRACTAdapter()))
    logger.info("Loaded ATtRACT adapter")
except Exception as e:
    logger.warning(f"Could not load ATtRACT adapter: {e}")

# --- circAtlas (Circular RNAs) ---
try:
    from template_package.adapters.circatlas_adapter import CircAtlasAdapter
    adapters.append(("circAtlas", CircAtlasAdapter()))
    logger.info("Loaded circAtlas adapter")
except Exception as e:
    logger.warning(f"Could not load circAtlas adapter: {e}")

# --- miRTarBase (Validated miRNA-Target Interactions) ---
try:
    from template_package.adapters.mirtarbase_adapter import MiRTarBaseAdapter
    adapters.append(("miRTarBase", MiRTarBaseAdapter()))
    logger.info("Loaded miRTarBase adapter")
except Exception as e:
    logger.warning(f"Could not load miRTarBase adapter: {e}")

# --- MatrixDB (Extracellular Matrix Interactions) ---
try:
    from template_package.adapters.matrixdb_adapter import MatrixDBAdapter
    adapters.append(("MatrixDB", MatrixDBAdapter()))
    logger.info("Loaded MatrixDB adapter")
except Exception as e:
    logger.warning(f"Could not load MatrixDB adapter: {e}")

# --- NLSdb / Signal Peptides ---
try:
    from template_package.adapters.nlsdb_adapter import NLSdbAdapter
    adapters.append(("NLSdb", NLSdbAdapter()))
    logger.info("Loaded NLSdb adapter")
except Exception as e:
    logger.warning(f"Could not load NLSdb adapter: {e}")

# --- NONCODE (Long Non-Coding RNAs) ---
try:
    from template_package.adapters.noncode_adapter import NONCODEAdapter
    adapters.append(("NONCODE", NONCODEAdapter()))
    logger.info("Loaded NONCODE adapter")
except Exception as e:
    logger.warning(f"Could not load NONCODE adapter: {e}")

# --- PsychENCODE GRN (Brain Gene Regulatory Network) ---
try:
    from template_package.adapters.psychencode_grn_adapter import PsychENCODEGRNAdapter
    adapters.append(("PsychENCODE_GRN", PsychENCODEGRNAdapter()))
    logger.info("Loaded PsychENCODE GRN adapter")
except Exception as e:
    logger.warning(f"Could not load PsychENCODE GRN adapter: {e}")

# --- MitoCarta (Mitochondrial Proteome) ---
try:
    from template_package.adapters.mitocarta_adapter import MitoCartaAdapter
    adapters.append(("MitoCarta", MitoCartaAdapter()))
    logger.info("Loaded MitoCarta adapter")
except Exception as e:
    logger.warning(f"Could not load MitoCarta adapter: {e}")

# --- UniProt Swiss-Prot (Human Proteome) ---
try:
    from template_package.adapters.uniprot_adapter import UniProtAdapter
    adapters.append(("UniProt", UniProtAdapter()))
    logger.info("Loaded UniProt adapter")
except Exception as e:
    logger.warning(f"Could not load UniProt adapter: {e}")

# --- GOA (Gene Ontology Annotations) ---
try:
    from template_package.adapters.goa_adapter import GOAAdapter
    adapters.append(("GOA", GOAAdapter()))
    logger.info("Loaded GOA adapter")
except Exception as e:
    logger.warning(f"Could not load GOA adapter: {e}")

# --- ClinVar (Clinical Variants) ---
try:
    from template_package.adapters.clinvar_adapter import ClinVarAdapter
    adapters.append(("ClinVar", ClinVarAdapter()))
    logger.info("Loaded ClinVar adapter")
except Exception as e:
    logger.warning(f"Could not load ClinVar adapter: {e}")

# --- REDIportal (RNA Editing Sites) ---
try:
    from template_package.adapters.rediportal_adapter import REDIportalAdapter
    adapters.append(("REDIportal", REDIportalAdapter()))
    logger.info("Loaded REDIportal adapter")
except Exception as e:
    logger.warning(f"Could not load REDIportal adapter: {e}")

# --- NeuroMorpho (Neuron Morphology) ---
try:
    from template_package.adapters.neuromorpho_adapter import NeuroMorphoAdapter
    adapters.append(("NeuroMorpho", NeuroMorphoAdapter()))
    logger.info("Loaded NeuroMorpho adapter")
except Exception as e:
    logger.warning(f"Could not load NeuroMorpho adapter: {e}")

# --- OmniPath (Signaling Network) ---
try:
    from template_package.adapters.omnipath_adapter import OmniPathAdapter
    adapters.append(("OmniPath", OmniPathAdapter()))
    logger.info("Loaded OmniPath adapter")
except Exception as e:
    logger.warning(f"Could not load OmniPath adapter: {e}")

# --- TISSUES (Tissue Expression) ---
try:
    from template_package.adapters.tissues_adapter import TISSUESAdapter
    adapters.append(("TISSUES", TISSUESAdapter()))
    logger.info("Loaded TISSUES adapter")
except Exception as e:
    logger.warning(f"Could not load TISSUES adapter: {e}")

# --- Bgee (Gene Expression) ---
try:
    from template_package.adapters.bgee_adapter import BgeeAdapter
    adapters.append(("Bgee", BgeeAdapter()))
    logger.info("Loaded Bgee adapter")
except Exception as e:
    logger.warning(f"Could not load Bgee adapter: {e}")

# ============================================================
# 3. Write nodes and edges from all adapters
# ============================================================
logger.info(f"Processing {len(adapters)} adapters...")

for name, adapter in adapters:
    try:
        # Collect nodes into a list to check if empty (BioCypher crashes on empty generators)
        logger.info(f"--- Writing nodes from {name} ---")
        nodes = list(adapter.get_nodes())
        if nodes:
            driver.write_nodes(iter(nodes))
        else:
            logger.info(f"--- {name} produced no nodes ---")

        logger.info(f"--- Writing edges from {name} ---")
        edges = list(adapter.get_edges())
        if edges:
            driver.write_edges(iter(edges))
        else:
            logger.info(f"--- {name} produced no edges ---")
        logger.info(f"--- Completed {name} ---")
    except Exception as e:
        logger.error(f"Error processing {name}: {e}")
        logger.warning(f"Skipping {name} due to error, continuing with remaining adapters...")

# ============================================================
# 4. Generate import call script
# ============================================================
driver.write_import_call()

# ============================================================
# 5. Summary
# ============================================================
logger.info("=" * 60)
logger.info("Import complete!")
logger.info(f"Processed {len(adapters)} adapters: {', '.join(n for n, _ in adapters)}")
logger.info("Check the 'biocypher-out' directory for CSVs and import scripts.")
logger.info("=" * 60)
