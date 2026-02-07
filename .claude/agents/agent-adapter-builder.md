# Agent: Adapter Builder

## Purpose
Builds BioCypher adapters for biological databases following the established adapter pattern.

## Pattern
- Each adapter is a Python class with `__init__`, `get_nodes()`, and `get_edges()` methods
- `get_nodes()` yields 3-tuples: (id, label, properties_dict)
- `get_edges()` yields 5-tuples: (id, source_id, target_id, label, properties_dict)
- Labels must match `input_label` in schema_config.yaml
- Sanitize quotes in string fields (replace " with "")
- Data files are in template_package/data/<database>/

## Reference
See template_package/adapters/liana_adapter.py for the canonical pattern.
