"""
Poison exon / NMD event Adapter for BioCypher.
Parses NMD (nonsense-mediated decay) events from supplementary XLSX data
(Nat Commun 2020, MOESM5). Each row links a gene to a poison-exon NMD event
at specific genomic coordinates.

Yields:
  - NMD event NODES (unique by gene + coordinates + AS type)
  - Gene-to-NMD-event EDGES (label: gene_has_nmd_event)
"""

import hashlib
import re
from pathlib import Path
from biocypher._logger import logger


# Pattern to detect Excel date-mangled gene names (e.g. SEPT2 -> 2019-09-02)
_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}")


class PoisonExonsAdapter:
    def __init__(self, data_dir="template_package/data/poison_exons"):
        self.data_dir = Path(data_dir)
        self.nmd_entries = []  # list of dicts with Gene, Coordinates, AS_type, Orphanet
        self._load_data()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _sanitize(self, text):
        if text is None:
            return ""
        text = str(text)
        text = text.replace('"', '""')
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        return text.strip()

    @staticmethod
    def _make_id(*parts):
        """Deterministic short hash from arbitrary string parts."""
        raw = "|".join(str(p) for p in parts)
        return hashlib.md5(raw.encode()).hexdigest()[:12]

    @staticmethod
    def _is_valid_gene(name):
        """Reject date-mangled gene names produced by Excel auto-formatting."""
        if not name:
            return False
        return not _DATE_RE.match(str(name))

    # ------------------------------------------------------------------
    # XLSX loading - openpyxl with zip/XML fallback
    # ------------------------------------------------------------------
    def _load_data(self):
        if not self.data_dir.exists():
            logger.warning("PoisonExonsAdapter: data directory not found")
            return

        main_file = self.data_dir / "41467_2020_17093_MOESM5_ESM.xlsx"
        if not main_file.exists():
            logger.warning("PoisonExonsAdapter: main XLSX not found")
            return

        try:
            self._load_with_openpyxl(main_file)
        except Exception as e1:
            logger.info(f"PoisonExonsAdapter: openpyxl failed ({e1}), trying XML fallback")
            try:
                self._load_with_xml_fallback(main_file)
            except Exception as e2:
                logger.warning(f"PoisonExonsAdapter: XML fallback also failed: {e2}")

        logger.info(f"PoisonExonsAdapter: Loaded {len(self.nmd_entries)} NMD events")

    def _load_with_openpyxl(self, path):
        """Primary parser using openpyxl."""
        import openpyxl
        wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
        ws = wb["NMD events"]

        header = None
        skipped = 0
        for row in ws.iter_rows(values_only=True):
            vals = list(row)
            # Detect the real header row (first cell == 'Gene')
            if header is None:
                if vals and str(vals[0]).strip() == "Gene":
                    header = [str(v).strip() if v else "" for v in vals]
                continue
            if not vals or vals[0] is None:
                continue
            rec = dict(zip(header, vals))
            gene = str(rec.get("Gene", "")).strip()
            if not self._is_valid_gene(gene):
                skipped += 1
                continue
            rec["Gene"] = gene
            self.nmd_entries.append(rec)
        wb.close()
        if skipped:
            logger.info(f"PoisonExonsAdapter: Skipped {skipped} date-mangled gene rows")

    def _load_with_xml_fallback(self, path):
        """Fallback: read the XLSX as a ZIP, parse shared strings + sheet XML."""
        import zipfile
        import xml.etree.ElementTree as ET

        z = zipfile.ZipFile(str(path))

        # 1. shared strings
        shared = []
        try:
            ss_xml = ET.parse(z.open("xl/sharedStrings.xml"))
            ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            for si in ss_xml.findall(".//s:si", ns):
                parts = []
                for t in si.iter():
                    if t.text:
                        parts.append(t.text)
                shared.append("".join(parts))
        except KeyError:
            pass

        # 2. find sheet1 (NMD events)
        ns_sheet = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        sheet_xml = ET.parse(z.open("xl/worksheets/sheet1.xml"))

        header = None
        skipped = 0
        for sheet_row in sheet_xml.findall(".//s:row", ns_sheet):
            cells = sheet_row.findall("s:c", ns_sheet)
            vals = []
            for c in cells:
                t_attr = c.get("t", "")
                v_el = c.find("s:v", ns_sheet)
                if v_el is not None and v_el.text is not None:
                    if t_attr == "s":
                        idx = int(v_el.text)
                        vals.append(shared[idx] if idx < len(shared) else "")
                    else:
                        vals.append(v_el.text)
                else:
                    vals.append(None)

            if header is None:
                if vals and str(vals[0]).strip() == "Gene":
                    header = [str(v).strip() if v else "" for v in vals]
                continue
            if not vals or vals[0] is None:
                continue
            # pad to header length
            while len(vals) < len(header):
                vals.append(None)
            rec = dict(zip(header, vals))
            gene = str(rec.get("Gene", "")).strip()
            if not self._is_valid_gene(gene):
                skipped += 1
                continue
            rec["Gene"] = gene
            self.nmd_entries.append(rec)
        z.close()
        if skipped:
            logger.info(f"PoisonExonsAdapter: Skipped {skipped} date-mangled gene rows (XML)")

    # ------------------------------------------------------------------
    # Parse coordinates helper
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_coords(coord_str):
        """Parse 'chr19:58353321-58353474:+' into components."""
        if not coord_str:
            return {}
        m = re.match(r"(chr\w+):(\d+)-(\d+):([+\-.])", str(coord_str))
        if m:
            return {
                "chromosome": m.group(1),
                "start": int(m.group(2)),
                "end": int(m.group(3)),
                "strand": m.group(4),
            }
        return {"raw_coordinates": str(coord_str)}

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self):
        """Yield NMD event nodes (unique by gene+coordinates+AS_type)."""
        seen = set()
        count = 0
        for rec in self.nmd_entries:
            gene = self._sanitize(rec.get("Gene"))
            coords = self._sanitize(rec.get("Coordinates"))
            as_type = self._sanitize(rec.get("AS type", ""))
            orphanet = self._sanitize(rec.get("Orphanet", ""))

            if not gene or not coords:
                continue

            node_id = f"nmd_event:{self._make_id(gene, coords, as_type)}"
            if node_id in seen:
                continue
            seen.add(node_id)

            coord_props = self._parse_coords(coords)
            props = {
                "gene": gene,
                "coordinates": coords,
                "as_type": as_type,
                "orphanet": orphanet,
                "source": "NatCommun_2020_poison_exons",
            }
            props.update(coord_props)

            count += 1
            yield node_id, "nmd_event", props

        logger.info(f"PoisonExonsAdapter: Yielded {count} NMD event nodes")

    def get_edges(self):
        """
        Yield edges: gene --[gene_has_nmd_event]--> nmd_event node.
        Source node = gene symbol, target = deterministic NMD-event id.
        """
        seen = set()
        count = 0
        for rec in self.nmd_entries:
            gene = self._sanitize(rec.get("Gene"))
            coords = self._sanitize(rec.get("Coordinates"))
            as_type = self._sanitize(rec.get("AS type", ""))
            orphanet = self._sanitize(rec.get("Orphanet", ""))

            if not gene or not coords:
                continue

            target_id = f"nmd_event:{self._make_id(gene, coords, as_type)}"
            edge_key = (gene, target_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)

            coord_props = self._parse_coords(coords)
            props = {
                "coordinates": coords,
                "as_type": as_type,
                "orphanet": orphanet,
                "source": "NatCommun_2020_poison_exons",
            }
            props.update(coord_props)

            count += 1
            yield (
                None,                       # auto edge id
                gene,                       # source: gene symbol
                target_id,                  # target: nmd event node
                "gene_has_nmd_event",       # label
                props,
            )

        logger.info(f"PoisonExonsAdapter: Yielded {count} edges")
