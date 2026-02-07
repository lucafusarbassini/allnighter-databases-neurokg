"""
eMouse Atlas (EMAPA) Adapter for BioCypher.

Parses the EMAPA ontology in OBO format and generates:
- AnatomicalTerm nodes (from [Term] stanzas)
- IsA relationship edges (from is_a lines)
- PartOf relationship edges (from relationship: part_of lines)

Data files:
- emapa.obo (EMAPA developmental anatomy ontology)
"""

import re
from pathlib import Path
from biocypher._logger import logger


class EMouseAdapter:
    def __init__(self, data_dir="template_package/data/emouse"):
        self.data_dir = Path(data_dir)
        self.terms = []
        self.is_a_edges = []
        self.part_of_edges = []
        self.other_relationships = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Parse all .obo files in the data directory."""
        obo_path = self.data_dir / 'emapa.obo'
        if not obo_path.exists():
            logger.warning("eMouse: emapa.obo not found")
            return
        self._parse_obo(obo_path)
        logger.info(
            f"eMouse: Loaded {len(self.terms)} terms, "
            f"{len(self.is_a_edges)} is_a edges, "
            f"{len(self.part_of_edges)} part_of edges, "
            f"{len(self.other_relationships)} other relationship edges"
        )

    def _parse_obo(self, fpath):
        """
        Parse OBO format file. Extracts:
        - id, name, namespace, def, synonyms, alt_ids, xrefs, comments
        - is_a relationships
        - part_of and other typed relationships
        - is_obsolete flag (skips obsolete terms)
        """
        current = None
        in_term = False

        with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.rstrip('\n')

                if line == '[Term]':
                    # Save previous term
                    if current and current.get('id') and not current.get('is_obsolete'):
                        self.terms.append(current)
                    current = {
                        'id': '',
                        'name': '',
                        'namespace': '',
                        'definition': '',
                        'synonyms': [],
                        'alt_ids': [],
                        'xrefs': [],
                        'comment': '',
                        'is_obsolete': False,
                    }
                    in_term = True
                    continue

                if line.startswith('[') and line.endswith(']'):
                    # Another stanza type (e.g. [Typedef]) -- save current term
                    if current and current.get('id') and not current.get('is_obsolete'):
                        self.terms.append(current)
                    current = None
                    in_term = False
                    continue

                if not in_term or current is None:
                    continue

                if line.startswith('id: '):
                    current['id'] = line[4:].strip()

                elif line.startswith('name: '):
                    current['name'] = line[6:].strip()

                elif line.startswith('namespace: '):
                    current['namespace'] = line[11:].strip()

                elif line.startswith('def: '):
                    # def: "definition text" [xrefs]
                    match = re.match(r'def:\s*"(.*?)"\s*\[', line)
                    if match:
                        current['definition'] = match.group(1)

                elif line.startswith('synonym: '):
                    # synonym: "text" SCOPE [xrefs]
                    match = re.match(r'synonym:\s*"(.*?)"\s+(\S+)', line)
                    if match:
                        current['synonyms'].append(match.group(1))

                elif line.startswith('alt_id: '):
                    current['alt_ids'].append(line[8:].strip())

                elif line.startswith('xref: '):
                    current['xrefs'].append(line[6:].strip())

                elif line.startswith('comment: '):
                    current['comment'] = line[9:].strip()

                elif line.startswith('is_obsolete: true'):
                    current['is_obsolete'] = True

                elif line.startswith('is_a: '):
                    # is_a: EMAPA:31859 ! polar body
                    parent_id = line[6:].split('!')[0].strip()
                    if parent_id and current['id']:
                        self.is_a_edges.append({
                            'child': current['id'],
                            'parent': parent_id,
                        })

                elif line.startswith('relationship: '):
                    # relationship: part_of EMAPA:36041 ! 1-cell stage conceptus
                    parts = line[14:].strip().split(None, 2)
                    if len(parts) >= 2:
                        rel_type = parts[0]
                        target_id = parts[1].split('!')[0].strip()
                        if target_id and current['id']:
                            if rel_type == 'part_of':
                                self.part_of_edges.append({
                                    'child': current['id'],
                                    'parent': target_id,
                                    'relationship': rel_type,
                                })
                            else:
                                self.other_relationships.append({
                                    'source': current['id'],
                                    'target': target_id,
                                    'relationship': rel_type,
                                })

        # Don't forget the last term
        if current and current.get('id') and not current.get('is_obsolete'):
            self.terms.append(current)

    def get_nodes(self):
        """
        Yield AnatomicalTerm nodes from EMAPA ontology.
        Each: (id, label, properties)
        """
        logger.info("eMouse: Generating AnatomicalTerm nodes...")
        count = 0

        for term in self.terms:
            term_id = term['id']
            props = {
                'name': self._sanitize(term['name']),
                'namespace': self._sanitize(term['namespace']),
                'definition': self._sanitize(term['definition'][:500] if term['definition'] else ''),
                'synonyms': '|'.join(self._sanitize(s) for s in term['synonyms']),
                'alt_ids': '|'.join(term['alt_ids']),
                'xrefs': '|'.join(term['xrefs']),
                'comment': self._sanitize(term['comment'][:300] if term['comment'] else ''),
                'source': 'EMAPA',
            }
            yield (f"emapa:{term_id}", "AnatomicalTerm", props)
            count += 1

        logger.info(f"eMouse: Generated {count} AnatomicalTerm nodes")

    def get_edges(self):
        """
        Yield is_a and part_of relationship edges.
        Each: (id, source, target, label, properties)
        """
        logger.info("eMouse: Generating relationship edges...")

        # is_a edges
        isa_count = 0
        for edge in self.is_a_edges:
            yield (
                f"emapa:isa:{edge['child']}_{edge['parent']}",
                f"emapa:{edge['child']}",
                f"emapa:{edge['parent']}",
                "IsA",
                {'relationship_type': 'is_a', 'source': 'EMAPA'},
            )
            isa_count += 1

        # part_of edges
        partof_count = 0
        for edge in self.part_of_edges:
            yield (
                f"emapa:partof:{edge['child']}_{edge['parent']}",
                f"emapa:{edge['child']}",
                f"emapa:{edge['parent']}",
                "PartOf",
                {'relationship_type': 'part_of', 'source': 'EMAPA'},
            )
            partof_count += 1

        # Other relationships (ends_at, starts_at, etc.)
        other_count = 0
        for edge in self.other_relationships:
            rel_type = edge['relationship']
            yield (
                f"emapa:{rel_type}:{edge['source']}_{edge['target']}",
                f"emapa:{edge['source']}",
                f"emapa:{edge['target']}",
                "OntologyRelationship",
                {'relationship_type': rel_type, 'source': 'EMAPA'},
            )
            other_count += 1

        logger.info(
            f"eMouse: Generated {isa_count} IsA + {partof_count} PartOf + "
            f"{other_count} other = {isa_count + partof_count + other_count} total edges"
        )
