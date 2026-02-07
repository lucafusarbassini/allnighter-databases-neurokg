"""
NeuroMorpho Adapter for BioCypher.

Loads NeuroMorpho.org human neuron morphology data and generates:
- NeuronMorphology nodes (digitally reconstructed neuron morphologies)

NeuroMorpho.org is the largest collection of publicly accessible
3D neuronal reconstructions, with detailed metadata on cell type,
brain region, species, and morphometric measurements.
"""

import json
from pathlib import Path
from biocypher._logger import logger


class NeuroMorphoAdapter:
    def __init__(self, data_dir="template_package/data/neuromorpho"):
        self.data_dir = Path(data_dir)
        self.neurons = []
        self._load_data()

    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    def _load_data(self):
        """Load NeuroMorpho human neuron data from JSON files."""
        json_files = sorted(self.data_dir.glob('neuromorpho_human*.json'))
        if not json_files:
            logger.warning("NeuroMorpho: No data files found")
            return

        logger.info(f"NeuroMorpho: Loading from {len(json_files)} JSON files...")
        count = 0
        seen = set()

        for jf in json_files:
            try:
                with open(jf, 'r', encoding='utf-8') as f:
                    data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"NeuroMorpho: Error reading {jf.name}: {e}")
                continue

            # Handle HAL-style API response
            if isinstance(data, dict) and '_embedded' in data:
                neurons = data['_embedded'].get('neuronResources', [])
            elif isinstance(data, list):
                neurons = data
            else:
                continue

            for neuron in neurons:
                nid = neuron.get('neuron_id')
                if nid is None or nid in seen:
                    continue
                seen.add(nid)

                brain_region = neuron.get('brain_region', [])
                if isinstance(brain_region, list):
                    brain_region = '; '.join(brain_region)

                cell_type = neuron.get('cell_type', [])
                if isinstance(cell_type, list):
                    cell_type = '; '.join(cell_type)

                self.neurons.append({
                    'neuron_id': str(nid),
                    'neuron_name': neuron.get('neuron_name', ''),
                    'archive': neuron.get('archive', ''),
                    'brain_region': brain_region,
                    'cell_type': cell_type,
                    'gender': neuron.get('gender', ''),
                    'age_classification': neuron.get('age_classification', ''),
                    'stain': neuron.get('stain', ''),
                    'protocol': neuron.get('protocol', ''),
                    'reconstruction_software': neuron.get('reconstruction_software', ''),
                    'soma_surface': neuron.get('soma_surface', ''),
                    'surface': neuron.get('surface', ''),
                    'volume': neuron.get('volume', ''),
                    'physical_integrity': neuron.get('physical_Integrity', ''),
                })
                count += 1

        logger.info(f"NeuroMorpho: Loaded {count} human neuron morphologies")

    def get_nodes(self):
        """
        Generate NeuronMorphology nodes.
        Yields: (id, label, properties)
        """
        logger.info("NeuroMorpho: Generating nodes...")
        count = 0

        for neuron in self.neurons:
            props = {
                'neuron_name': self._sanitize(neuron['neuron_name']),
                'archive': self._sanitize(neuron['archive']),
                'brain_region': self._sanitize(neuron['brain_region']),
                'cell_type': self._sanitize(neuron['cell_type']),
                'gender': self._sanitize(neuron['gender']),
                'age_classification': self._sanitize(neuron['age_classification']),
                'stain': self._sanitize(neuron['stain']),
                'protocol': self._sanitize(neuron['protocol']),
                'reconstruction_software': self._sanitize(neuron['reconstruction_software']),
                'soma_surface': self._sanitize(neuron['soma_surface']),
                'surface': self._sanitize(neuron['surface']),
                'volume': self._sanitize(neuron['volume']),
                'physical_integrity': self._sanitize(neuron['physical_integrity']),
                'source': 'NeuroMorpho',
            }

            yield (f"neuromorpho:{neuron['neuron_id']}", "NeuronMorphology", props)
            count += 1

        logger.info(f"NeuroMorpho: Generated {count} NeuronMorphology nodes")

    def get_edges(self):
        """No edges."""
        logger.info("NeuroMorpho: No edges to generate")
        return iter([])
