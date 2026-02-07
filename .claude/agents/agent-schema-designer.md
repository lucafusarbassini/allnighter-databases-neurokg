# Agent: Schema Designer

## Purpose
Designs and updates the BioCypher schema configuration for new database integrations.

## Guidelines
- Use BioLink ontology names in sentence case for schema keys
- Use PascalCase for input_label values
- Use `is_a` for custom ontology branches
- Properties: str, int, float, bool, str[], int[], float[], bool[]
- Cross-references stored as JSON string
- Reference: config/schema_config.yaml and INSTRUCTIONS.md
