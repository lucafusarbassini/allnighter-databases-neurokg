# TODO - Current Work Items

## Completed
1. [DONE] All 12 remaining data directories now have adapters:
   - 3DGenome, Adhesome, BRENDA, ChapNet, ENCORI, EpiMap
   - EVpedia, LipidOntology (LION), RaftProt, SILVA, TarBase/miRDB, UniCarbKB
2. [DONE] Schema config updated (102 entries)
3. [DONE] create_knowledge_graph.py updated (65 adapters)

## Next Steps
1. [ ] Run end-to-end pipeline test (docker compose up)
2. [ ] Verify all adapters load successfully
3. [ ] Check for cross-reference linking opportunities between databases
4. [ ] Add more neuroscience databases from task.md list:
   - mousebrain.org (single-cell atlas)
   - Allen Brain Atlas ABC Atlas (additional data)
   - eMouse atlas
   - SynGO portal (additional data)
   - PsychENCODE (additional data layers)
5. [ ] Consider adding adapters for databases referenced by papers/URLs in task.md
   that don't yet have downloaded data
6. [ ] Optimize adapter loading for large datasets (EpiMap, TarBase)
7. [ ] Add more cross-linking edges between entities from different databases
