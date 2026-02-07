"""
Allen Brain Atlas Adapter for BioCypher.

Loads Allen Brain Atlas structure graph (mouse brain anatomy hierarchy)
and generates:
- BrainRegion nodes (anatomical brain structures)
- BrainRegionPartOf edges (parent-child hierarchy)
"""

import json
import re
from pathlib import Path
from biocypher._logger import logger


class AllenBrainAdapter:
    def __init__(self, data_dir="template_package/data/allen_brain"):
        self.data_dir = Path(data_dir)
        self.structures = {}
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = re.sub(r'<[^>]+>', '', text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load Allen Brain structure graph."""
        sg_path = self.data_dir / 'structure_graph.json'
        if sg_path.exists() and sg_path.stat().st_size > 100:
            logger.info("AllenBrain: Loading structure graph...")
            with open(sg_path, 'r') as f:
                data = json.load(f)

            if data.get('success') and data.get('msg'):
                # The msg is a list with one root element containing nested children
                for root in data['msg']:
                    self._flatten_structure(root, parent_id=None)

        logger.info(f"AllenBrain: Loaded {len(self.structures)} brain structures")

    def _flatten_structure(self, node, parent_id):
        """Recursively flatten the nested structure tree."""
        struct_id = node.get('id')
        if struct_id is None:
            return

        self.structures[struct_id] = {
            'id': struct_id,
            'acronym': node.get('acronym', ''),
            'name': node.get('name', ''),
            'color_hex': node.get('color_hex_triplet', ''),
            'graph_order': node.get('graph_order', 0),
            'st_level': node.get('st_level', 0),
            'parent_id': parent_id,
            'atlas_id': node.get('atlas_id'),
            'hemisphere_id': node.get('hemisphere_id', 3),
        }

        children = node.get('children', [])
        if children:
            for child in children:
                self._flatten_structure(child, parent_id=struct_id)

    def get_nodes(self):
        """
        Generate brain region nodes.
        Yields: (id, label, properties)
        """
        logger.info("AllenBrain: Generating nodes...")
        count = 0

        for struct_id, struct in self.structures.items():
            node_id = f"ABA:{struct_id}"
            props = {
                'name': self._sanitize(struct.get('name', '')),
                'acronym': self._sanitize(struct.get('acronym', '')),
                'color_hex': struct.get('color_hex', ''),
                'graph_order': struct.get('graph_order', 0),
                'structure_level': struct.get('st_level', 0),
                'source': 'AllenBrainAtlas',
            }

            yield (node_id, "BrainRegion", props)
            count += 1

        logger.info(f"AllenBrain: Generated {count} BrainRegion nodes")

    def get_edges(self):
        """
        Generate brain region hierarchy edges.
        Yields: (id, source, target, label, properties)
        """
        logger.info("AllenBrain: Generating edges...")
        count = 0

        for struct_id, struct in self.structures.items():
            parent_id = struct.get('parent_id')
            if parent_id is not None and parent_id in self.structures:
                yield (
                    None,
                    f"ABA:{struct_id}",
                    f"ABA:{parent_id}",
                    "BrainRegionPartOf",
                    {}
                )
                count += 1

        logger.info(f"AllenBrain: Generated {count} BrainRegionPartOf edges")
