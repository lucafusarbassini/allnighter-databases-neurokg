"""
Reactome Pathway Adapter for BioCypher.

Loads Reactome pathway data and generates:
- ReactomePathway nodes
- GeneParticipatesInPathway edges (Gene → ReactomePathway)
- PathwayIsPartOfPathway edges (pathway hierarchy)
"""

from pathlib import Path
from biocypher._logger import logger


class ReactomeAdapter:
    def __init__(self, data_dir="template_package/data/reactome"):
        self.data_dir = data_dir
        self.pathways = {}  # R-HSA-XXXX -> {name, url}
        self.gene_pathway_links = []  # (uniprot_id, pathway_id, evidence)
        self.pathway_hierarchy = []  # (parent_id, child_id)
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Reactome UniProt-to-pathway mapping (human only)."""
        filepath = Path(self.data_dir) / 'UniProt2Reactome_All_Levels.txt'
        if not filepath.exists():
            logger.warning("Reactome: UniProt2Reactome file not found")
            return

        logger.info("Reactome: Loading UniProt-to-pathway mappings (human only)...")
        total = 0
        kept = 0

        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                total += 1
                uniprot_id = parts[0].strip()
                pathway_id = parts[1].strip()
                url = parts[2].strip()
                pathway_name = self._sanitize(parts[3].strip())
                evidence = parts[4].strip()
                species = parts[5].strip()

                # Filter: only human pathways (R-HSA-*)
                if species != 'Homo sapiens':
                    continue
                if not pathway_id.startswith('R-HSA-'):
                    continue

                kept += 1

                # Register pathway
                if pathway_id not in self.pathways:
                    self.pathways[pathway_id] = {
                        'name': pathway_name,
                        'url': url,
                    }

                # Register gene-pathway link
                self.gene_pathway_links.append({
                    'uniprot_id': uniprot_id,
                    'pathway_id': pathway_id,
                    'evidence': evidence,
                })

        logger.info(f"Reactome: Loaded {len(self.pathways)} human pathways, "
                     f"{kept} gene-pathway links (from {total} total)")

        # Load pathway hierarchy
        self._load_pathway_hierarchy()

    def _load_pathway_hierarchy(self):
        """Load ReactomePathwaysRelation.txt for pathway parent-child relationships."""
        rel_path = Path(self.data_dir) / 'ReactomePathwaysRelation.txt'
        if not rel_path.exists():
            logger.info("Reactome: No pathway hierarchy file found")
            return

        logger.info("Reactome: Loading pathway hierarchy...")
        count = 0

        with open(rel_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                parent_id = parts[0].strip()
                child_id = parts[1].strip()

                # Only human pathways
                if parent_id.startswith('R-HSA-') and child_id.startswith('R-HSA-'):
                    self.pathway_hierarchy.append((parent_id, child_id))
                    count += 1

                    # Ensure both pathways exist in our map
                    for pid in [parent_id, child_id]:
                        if pid not in self.pathways:
                            self.pathways[pid] = {
                                'name': '',
                                'url': f'https://reactome.org/content/detail/{pid}',
                            }

        # Enrich pathway names from ReactomePathways.txt if available
        names_path = Path(self.data_dir) / 'ReactomePathways.txt'
        if names_path.exists():
            with open(names_path, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        pid = parts[0].strip()
                        pname = self._sanitize(parts[1].strip())
                        species = parts[2].strip()
                        if pid.startswith('R-HSA-') and pid in self.pathways:
                            if not self.pathways[pid]['name']:
                                self.pathways[pid]['name'] = pname

        logger.info(f"Reactome: Loaded {count} pathway hierarchy relations")

    def get_nodes(self):
        """
        Generate ReactomePathway nodes.
        Yields: (id, label, properties)
        """
        logger.info("Reactome: Generating pathway nodes...")
        count = 0

        for pathway_id, data in self.pathways.items():
            props = {
                'name': data['name'],
                'url': data['url'],
                'source': 'Reactome',
            }
            yield (pathway_id, "ReactomePathway", props)
            count += 1

        logger.info(f"Reactome: Generated {count} ReactomePathway nodes")

    def get_edges(self):
        """
        Generate GeneParticipatesInPathway edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("Reactome: Generating gene-pathway edges...")
        count = 0
        seen = set()

        for link in self.gene_pathway_links:
            # Deduplicate: same gene-pathway pair may appear with different evidence
            key = (link['uniprot_id'], link['pathway_id'])
            if key in seen:
                continue
            seen.add(key)

            props = {
                'evidence': link['evidence'],
                'source': 'Reactome',
            }

            yield (
                None,
                link['uniprot_id'],
                link['pathway_id'],
                "GeneParticipatesInPathway",
                props
            )
            count += 1

        logger.info(f"Reactome: Generated {count} GeneParticipatesInPathway edges")

        # Pathway hierarchy edges
        hier_count = 0
        for parent_id, child_id in self.pathway_hierarchy:
            yield (
                None,
                child_id,
                parent_id,
                "PathwayIsPartOfPathway",
                {'source': 'Reactome'}
            )
            hier_count += 1

        logger.info(f"Reactome: Generated {hier_count} PathwayIsPartOfPathway edges")
