
The goal is to download and process a large set of biologically important databases related to cell biology and neuroscience. Specifically, this repo describes how to adapt databases, using the Liana database as an example, to the BioCypher and BioLink framework, respecting the BioLink ontology and allowing the construction of knowledge graphs. The goal of this night's work is therefore to download all the databases and harmonize them so they can then be combined into a single, enormous knowledge graph. One area where there is some arbitrariness to disambiguate is the ontology. In the case of Liana, for example, we saw that the database contained proteins represented either as stand-alone proteins or as quaternary structures, thus with multiple gene names. This must be inserted into the BioLink ontological schema. Once this is done, it is then accessible for searches in subsequent databases. For this reason, I suggest processing all databases sequentially, one after the other, so that the ever-moderate and controlled growth of the ontology can be used to analyze subsequent databases. Thus, a database is processed and may occasionally create further branches in the ontology, which the subsequent database must be able to see and use. Some of them you processed in previous iterations and are available here; your task is to complete the extraction to finally include ALL databases of my list.


These are the databases:

AREA 1: MOLECULAR AND CELLULAR BIOLOGY
CHEBI: https://www.ebi.ac.uk/chebi/
LIPIDMAPS: https://www.lipidmaps.org/
Rhea: https://www.rhea-db.org/help/download
Lipid Ontology: http://www.lipidontology.com/
KEGG: https://www.kegg.jp/kegg/pathway.html
https://academic.oup.com/nar/article/53/D1/D966/7905300
https://rnasysu.com/encori/index.php
https://academic.oup.com/nar/article/50/D1/D287/6353804
https://maayanlab.cloud/chea3/
https://www.proteinatlas.org/
https://academic.oup.com/nar/article/52/D1/D174/7420101
https://www.ebi.ac.uk/interpro/
https://www.nature.com/articles/s41467-024-52146-3
https://pubmed.ncbi.nlm.nih.gov/16397007/
https://pubmed.ncbi.nlm.nih.gov/34079125/
https://academic.oup.com/nar/article/49/D1/D1541/5974091
https://compartments.jensenlab.org/Search
https://www.cell.com/biophysj/fulltext/S0006-3495(24)03067-4
https://www.proteinatlas.org/humanproteome/subcellular (i love this one)
https://academic.oup.com/nar/article/48/D1/D174/5588346
https://ptmcode.embl.de/index.cgi
https://www.nature.com/articles/s41586-020-2188-x
https://www.nature.com/articles/s41467-018-08191-w
https://www.pnas.org/doi/10.1073/pnas.1916584117?utm_source=chatgpt.com
https://pmc.ncbi.nlm.nih.gov/articles/PMC11572202/?utm_source=chatgpt.com
https://pubmed.ncbi.nlm.nih.gov/38459036/
https://pmc.ncbi.nlm.nih.gov/articles/PMC11404449/?utm_source=chatgpt.com
https://alphafoldserver.com/welcome - placeholder, to say that maaaany binding events are probably available as predictions across several molecule classes now
https://pubmed.ncbi.nlm.nih.gov/29106644/
https://pubmed.ncbi.nlm.nih.gov/26339475/
https://academic.oup.com/nar/article/45/D1/D658/2333932?utm_source=chatgpt.com&login=true
https://academic.oup.com/bioinformatics/article/36/12/3941/5824292
https://doi.org/10.1093/nar/gkae1005
https://sabiork.h-its.org/
https://www.brenda-enzymes.org/index.php
https://research.bioinformatics.udel.edu/iptmnet/?utm_source=chatgpt.com
https://academic.oup.com/nar/article/53/D1/D377/7889255?login=true
https://screen.wenglab.org/
https://fantom.gsc.riken.jp/
https://www.nature.com/articles/s41594-022-00910-8?fromPaywallRec=false
https://www.ebi.ac.uk/complexportal/home
https://www.ebi.ac.uk/intact/home
https://academic.oup.com/nar/article/50/D1/D54/6424762?login=true
https://comppi.linkgroup.hu/?utm_source=chatgpt.com
https://gpcrdb.org/?utm_source=chatgpt.com
https://translocatome.linkgroup.hu/
https://pubmed.ncbi.nlm.nih.gov/22782549/
https://academic.oup.com/nar/article/52/D1/D1694/7416396
https://www.nature.com/articles/s41467-020-17093-9 contains a large list of poison exons and other unproductive splicing events (suppl data)
https://www.nature.com/articles/s41586-025-08878-3
http://cpdb.molgen.mpg.de/ old school but might have something useful
https://pmc.ncbi.nlm.nih.gov/articles/PMC3998127
https://academic.oup.com/nar/article/47/D1/D221/5160993
https://www.nature.com/articles/s41588-019-0500-1
https://3dgenome.fsm.northwestern.edu/?utm_source=chatgpt.com
https://academic.oup.com/nar/article/51/D1/D159/6754910
https://alphafold.ebi.ac.uk/ not really for now as 3d structure (and variable conformation) can surely be represented with graphs, but it is low priority and likely a project by itself
https://pmc.ncbi.nlm.nih.gov/articles/PMC11655564/?utm_source=chatgpt.com probably has good references for interesting databases, not urgent
http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index
https://www.nature.com/articles/s41586-020-2493-4
https://academic.oup.com/nar/article/50/D1/D222/6446528
https://pubmed.ncbi.nlm.nih.gov/33693667/
https://arrestindb.org/
https://academic.oup.com/nar/article/40/D1/D1241/2903429?utm_source=chatgpt.com&login=true
https://zhanggroup.org/
https://doi.org/10.1093/database/baw035
https://academic.oup.com/nar/article/50/D1/D231/6454677
https://www.tcdb.org/?utm_source=chatgpt.com
https://pathguide.org/ resource of other, secondary importance databases
https://www.expasy.org/ resource of other, secondary importance databases
https://academic.oup.com/nar/article/53/D1/D1677/7903369?utm_source=chatgpt.com&login=true
https://matrisomedb.org/
https://www.science.org/doi/10.1126/scisignal.aaz0274?utm_source=chatgpt.com
https://www.proteinatlas.org/humanproteome/tissue/secretome?utm_source=chatgpt.com
https://academic.oup.com/nar/article/43/D1/D1140/2437426
https://zhanggroup.org/BioLiP/download.html?utm_source=chatgpt.com
https://metalpdb.cerm.unifi.it/?utm_source=chatgpt.com
https://ngdc.cncb.ac.cn/databasecommons/database/id/5672?utm_source=chatgpt.com
https://academic.oup.com/nar/article/48/D1/D288/5613676?utm_source=chatgpt.com&login=true
https://db.phasep.pro/?utm_source=chatgpt.com
https://degronopedia.com/?utm_source=chatgpt.com
http://elm.eu.org/
https://mobidb.org/?utm_source=chatgpt.com
https://trep-db.uzh.ch/
https://www.dfam.org/
https://peptideatlas.org/repository/?utm_source=chatgpt.com
https://pubmed.ncbi.nlm.nih.gov/34986596/
https://openprot.org/
https://pubmed.ncbi.nlm.nih.gov/26527729/
https://compbio.mit.edu/epimap/
https://academic.oup.com/database/article/doi/10.1093/database/baae007/7606618?utm_source=chatgpt.com&login=true
https://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?g=encodeCcreCombined&token=0.tVseYbX55WGfkv-5luSofJnZ5Lu4VyVNDCOSnJv-BoL--7EhpmXaM8pkUNpLXfIpSrvCIGZ7Sp9NeGEf31UmnzL-tO6qwm1KxnUjhwDntMN_pzXnJevXUBPMs8pqdut434NRsaNMNnYsx31qCjA4NZLXOBJKk1WXrt6WUQf4rBGiTobxtiZ85l1qQrhSnNzTWmjHGhip4f8cx5YP7ji3dPdE86eZoDoMvV3-27lUMhaiL_1VTq-3nBfW4tyyDp4f9h0FhI8J8VlTf7rxrdmKcE_A8nTc1iAduCEKlPKyJABZwNnocbHXq_Eu5kc5A49qp1yhln69shG8Kgfjl7lNcHVeWJKx_Eq1PIC-APm99TdhE0DtQ2S81lffXi0B-pioTgL6l3FovNY2r6-xBCN72woZpJH1uw7Kksmw3LNwuB1T5dkxuwY7WTyHsLk39JwHSbmaO4mcznjuKE4-sZL6uaQoRw9FQkqmWcncBTDXsnU1s7fUdo6CDoySAQ8bh4awtzQEA1BKm4xhBRYLsuFppd-T005ILVCccZhv_fPMFT8wgHAd4qfpxuqH1w2obskEOgCLepfQT-YH2uo1CRdT6_HB1yinkjiTKfWVCqflU2j73EbGiyRdy5uJA3VEU3n_aNqeBIL4JelbfTbKnXwjr0eK0nvxmwB3RgvGwv4UY75NstgEnSjV_SDsf4r1IfGeqx5_gE0MeQWitJsicBePFrpqpQMKkCM_p0AlxWr4RW4TrIygMIKuDfVfd7OTyrH_Y8F_cjK9PcUr9B0cDVvu--5sjLlqmDETC-cEYDwMXGZCXCQJRPAbXp6B09LOxsZayPwMKAHk6S-VuNgLF4hvNiT6BPU46yD2GSJkImpGmYJJpg821wmI7Rl80J-iEMTDkJQU1m3zQY74oAsVqG-bt32ZCddQ6PEb8-W_2GLxojE.230bZvQyfeWmY2wcW-4-PA.c2083865027dad77d4b93e3e42ac3fba9ae91d99f36e803f5f479d0a46988f9d
https://membranome.org/
https://evpedia.info/evpedia2_xe/
http://www.exocarta.org/
https://www.ebi.ac.uk/pdbe/pisa/?utm_source=chatgpt.com
https://rnagranuledb.lunenfeld.ca/?utm_source=chatgpt.com
https://raftprot.org/?utm_source=chatgpt.com
https://www.arb-silva.de/?utm_source=chatgpt.com
https://gtrnadb.org/?utm_source=chatgpt.com
https://pubmed.ncbi.nlm.nih.gov/37889077/
https://academic.oup.com/nar/article/49/D1/D165/5983627?utm_source=chatgpt.com&login=true
https://academic.oup.com/nar/article/52/D1/D52/7280546?utm_source=chatgpt.com&login=true
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6541-0?utm_source=chatgpt.com
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4?utm_source=chatgpt.com
https://adhesome.org/?utm_source=chatgpt.com
https://netbio.bgu.ac.il/chapnet/?utm_source=chatgpt.com
https://slc.bioparadigms.org/?utm_source=chatgpt.com
http://www.zhounan.org/ferrdb/current/
https://repairtoire.genesilico.pl/?utm_source=chatgpt.com
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04239-9?utm_source=chatgpt.com
https://www.ncbi.nlm.nih.gov/research/histonedb/?utm_source=chatgpt.com
https://academic.oup.com/nar/article/45/D1/D750/2770646?utm_source=chatgpt.com&login=true
https://www.unicarbkb.org/
https://www.nature.com/articles/s41556-024-01435-6
https://academic.oup.com/database/article/doi/10.1093/database/baae040/7683720?utm_source=chatgpt.com&login=true
https://www.nature.com/articles/s41587-021-01006-2
https://www.sciencedirect.com/science/article/pii/S0896627314004486?via%3Dihub

AREA 2: NEUROSCIENCE
http://mousebrain.org/
https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas
https://connectivity.brain-map.org/
https://www.nature.com/articles/s41556-021-00787-7 general cell ontology (also https://github.com/obophenotype/cell-ontology )
Take NeuroKG and see what it contains yet
Allen anatomical hierarchy (to be provided by Luca)
eMouse atlas: https://www.emouseatlas.org/emap/home.html
http://bag.linked-brain-data.org/
https://beta.brainkb.org/ contact them?
https://www.psychencode.org/
https://www.syngoportal.org/
https://modeldb.science/
https://neuroelectro.org/index.html (maybe there’s smth newer)
[you should next extend with more relevant neuro databases, and add them here to the list!]
