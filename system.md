# System Discovery

## Environment
- OS: Linux 6.8.0-52-generic
- Python: Available (^3.10 required)
- Package manager: Poetry (pyproject.toml)
- Git: Available
- Docker: Available (docker-compose.yml present)

## Dependencies (current)
- biocypher==0.10.1
- pyarrow>=22.0.0
- python ^3.10

## Key Paths
- Project root: /workspace
- Adapters: /workspace/template_package/adapters/
- Data: /workspace/template_package/data/
- Config: /workspace/config/
- Scripts: /workspace/scripts/
- Main entry: /workspace/create_knowledge_graph.py
- Schema: /workspace/config/schema_config.yaml

## Available Tools
- Python with pandas, biocypher, pyarrow
- Git for version control
- Docker for deployment pipeline
- pip/poetry for dependency management
