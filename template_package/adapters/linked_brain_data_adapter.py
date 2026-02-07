"""
Linked Brain Data Adapter for BioCypher.

Loads neuroscience entity data from the Linked Brain Data project
(linked-neuron-data.org / BrainKnow) and generates:
- NamedThing nodes representing neuroscience entities (diseases, brain
  regions, neuron types, genes, transmitters, pathways)

The data comes in two formats:
1. brainknow_entities.json -- list of dicts with uri, name, category_color,
   count, and synonyms fields.
2. search_*.json -- list of lists [uri, name, color, count, synonym_text]
   representing search results for specific queries.

Entities are deduplicated by URI across all files.

Expected data directory: template_package/data/linked_brain_data/
Expected files: brainknow_entities.json, search_*.json
"""

import json
from pathlib import Path
from biocypher._logger import logger

# Map the color codes used in the Linked Brain Data UI to entity categories
_COLOR_CATEGORY = {
    "lavender": "disease",
    "lightpink": "brain region",
    "CornflowerBlue": "neuron",
    "cornflowerblue": "neuron",
    "Aquamarine": "gene",
    "aquamarine": "gene",
    "NavajoWhite": "pathway",
    "navajowhite": "pathway",
    "Salmon": "transmitter",
    "salmon": "transmitter",
}


class LinkedBrainDataAdapter:
    def __init__(self, data_dir="template_package/data/linked_brain_data"):
        self.data_dir = Path(data_dir)
        self.entities = {}  # uri -> dict, deduplicated
        self._load_data()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _sanitize(self, text):
        """Sanitize string for CSV safety."""
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
        return text.strip()

    @staticmethod
    def _infer_category_from_uri(uri):
        """Extract entity category from the URI path structure."""
        # URIs like http://www.linked-neuron-data.org/resource/<category>/<name>
        prefix = "http://www.linked-neuron-data.org/resource/"
        if uri.startswith(prefix):
            remainder = uri[len(prefix):]
            parts = remainder.split("/")
            if parts:
                return parts[0]
        return "unknown"

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def _load_data(self):
        """Load all JSON data files from the linked_brain_data directory."""
        if not self.data_dir.exists():
            logger.warning(
                f"LinkedBrainData: Data directory not found: {self.data_dir}"
            )
            return

        json_files = sorted(self.data_dir.glob("*.json"))
        if not json_files:
            logger.warning(
                f"LinkedBrainData: No JSON files found in {self.data_dir}"
            )
            return

        for fpath in json_files:
            try:
                with open(fpath, "r", encoding="utf-8") as f:
                    data = json.load(f)

                if not isinstance(data, list) or len(data) == 0:
                    continue

                # Detect format: list of dicts vs list of lists
                first = data[0]
                if isinstance(first, dict):
                    self._parse_dict_format(data, fpath.name)
                elif isinstance(first, list):
                    self._parse_list_format(data, fpath.name)
                else:
                    logger.warning(
                        f"LinkedBrainData: Unknown format in {fpath.name}"
                    )

            except Exception as e:
                logger.warning(
                    f"LinkedBrainData: Error reading {fpath.name}: {e}"
                )

        logger.info(
            f"LinkedBrainData: Loaded {len(self.entities)} unique entities"
        )

    def _parse_dict_format(self, data, filename):
        """Parse brainknow_entities.json format (list of dicts)."""
        for entry in data:
            uri = entry.get("uri", "")
            if not uri:
                continue
            if uri in self.entities:
                continue  # deduplicate

            color = entry.get("category_color", "")
            category = _COLOR_CATEGORY.get(color, "")
            if not category:
                category = self._infer_category_from_uri(uri)

            synonyms_raw = entry.get("synonyms", "")
            synonyms = ""
            if synonyms_raw and synonyms_raw.startswith("Synonym: "):
                synonyms = synonyms_raw[len("Synonym: "):]

            self.entities[uri] = {
                "uri": uri,
                "name": entry.get("name", ""),
                "category": category,
                "count": entry.get("count", 0),
                "synonyms": synonyms,
                "source_file": filename,
            }

    def _parse_list_format(self, data, filename):
        """Parse search_*.json format (list of lists)."""
        for entry in data:
            if not isinstance(entry, list) or len(entry) < 4:
                continue

            uri = entry[0] if len(entry) > 0 else ""
            name = entry[1] if len(entry) > 1 else ""
            color = entry[2] if len(entry) > 2 else ""
            count = entry[3] if len(entry) > 3 else 0
            synonym_text = entry[4] if len(entry) > 4 else ""

            if not uri:
                continue
            if uri in self.entities:
                continue  # deduplicate

            category = _COLOR_CATEGORY.get(color, "")
            if not category:
                category = self._infer_category_from_uri(uri)

            synonyms = ""
            if synonym_text and synonym_text.startswith("Synonym: "):
                synonyms = synonym_text[len("Synonym: "):]

            self.entities[uri] = {
                "uri": uri,
                "name": name,
                "category": category,
                "count": count,
                "synonyms": synonyms,
                "source_file": filename,
            }

    # ------------------------------------------------------------------
    # Node generator
    # ------------------------------------------------------------------
    def get_nodes(self):
        """
        Generate NamedThing nodes from Linked Brain Data entities.

        Each entity carries its category, reference count, and synonyms.

        Yields: (id, label, properties)
        """
        logger.info(
            f"LinkedBrainData: Generating nodes from "
            f"{len(self.entities)} entities..."
        )
        count = 0

        for uri, entity in self.entities.items():
            node_id = self._sanitize(uri)
            name = self._sanitize(entity["name"])
            category = self._sanitize(entity["category"])

            props = {
                "name": name,
                "category": category,
                "reference_count": entity["count"],
                "source": "LinkedBrainData",
            }
            if entity["synonyms"]:
                props["synonyms"] = self._sanitize(entity["synonyms"])

            yield (node_id, "named thing", props)
            count += 1

        logger.info(f"LinkedBrainData: Generated {count} NamedThing nodes")

    # ------------------------------------------------------------------
    # Edge generator
    # ------------------------------------------------------------------
    def get_edges(self):
        """
        Linked Brain Data search results do not encode explicit
        relationships. This method yields nothing; edges could be added
        in the future when relationship data is available.
        """
        logger.info(
            "LinkedBrainData: No edges (entity data only)"
        )
        return
        yield  # makes this a generator
