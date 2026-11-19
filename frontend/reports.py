"""Report generation utilities (PDF, Word DOCX, and ZIP archive) using pure standard-library Python."""
from __future__ import annotations

import datetime
import io
import zipfile
from xml.sax.saxutils import escape


def _pdf_escape(text: str) -> str:
    """Normalize scientific symbols and escape special characters for a PDF literal string."""
    normalized = (
        str(text)
        .replace("Å²", "A^2")
        .replace("Å", "A")
        .replace("—", "-")
        .replace("•", "*")
    )
    cleaned = normalized.encode("ascii", errors="replace").decode("ascii")
    return cleaned.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def generate_pdf_report(
    project_title: str,
    parameters: dict,
    metrics: list[dict],
    table_rows: list[dict],
) -> bytes:
    """Compile project parameters, metrics, and table rows into a valid standard PDF document."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: list[tuple[str, int, bool]] = [
        (f"{project_title} - Structural Biology & Docking Report", 15, True),
        (f"Generated: {timestamp}", 10, False),
        ("", 10, False),
        ("1. Project & Execution Parameters", 13, True),
    ]
    for key, val in parameters.items():
        lines.append((f"  - {key}: {val}", 10, False))

    lines.append(("", 10, False))
    lines.append(("2. RDKit Physicochemical & Lipinski Metrics", 13, True))
    for item in metrics:
        title = item.get("title", "Metric")
        value = item.get("value", "")
        status = item.get("status", "")
        lines.append((f"  - {title}: {value} [{status}]", 10, False))

    lines.append(("", 10, False))
    lines.append(("3. AutoDock Vina v1.2.7 Ranked Binding Modes", 13, True))
    if table_rows:
        headers = list(table_rows[0].keys())
        lines.append(("  " + " | ".join(str(h) for h in headers), 8, True))
        lines.append(("  " + "-" * 85, 8, False))
        for row in table_rows:
            row_str = " | ".join(str(row.get(h, "")) for h in headers)
            lines.append(("  " + row_str[:105], 8, False))
    else:
        lines.append(("  No table rows recorded.", 10, False))

    content_ops = ["BT", "40 790 Td"]
    y_pos = 790
    for text, size, is_bold in lines:
        if y_pos < 50:
            break
        font_key = "/F2" if is_bold else "/F1"
        step = size + 5
        content_ops.append(f"{font_key} {size} Tf")
        content_ops.append(f"({_pdf_escape(text)}) Tj")
        content_ops.append(f"0 -{step} Td")
        y_pos -= step
    content_ops.append("ET")
    stream_bytes = "\n".join(content_ops).encode("ascii")

    objects = [
        b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n",
        b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
        (
            b"3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 842] "
            b"/Resources << /Font << /F1 4 0 R /F2 5 0 R >> >> /Contents 6 0 R >>\nendobj\n"
        ),
        b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
        b"5 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>\nendobj\n",
        f"6 0 obj\n<< /Length {len(stream_bytes)} >>\nstream\n".encode("ascii")
        + stream_bytes
        + b"\nendstream\nendobj\n",
    ]

    pdf = bytearray(b"%PDF-1.4\n")
    offsets: list[int] = []
    for obj in objects:
        offsets.append(len(pdf))
        pdf.extend(obj)

    xref_pos = len(pdf)
    pdf.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    pdf.extend(b"0000000000 65535 f \n")
    for off in offsets:
        pdf.extend(f"{off:010d} 00000 n \n".encode("ascii"))

    pdf.extend(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode(
            "ascii"
        )
    )
    return bytes(pdf)


def _docx_paragraph(text: str, bold: bool = False, color: str = "0A0F1C", size_half_pt: int = 22) -> str:
    """Return an OpenXML paragraph element."""
    b_tag = "<w:b/>" if bold else ""
    return (
        "<w:p>"
        f"<w:r><w:rPr>{b_tag}<w:color w:val=\"{color}\"/><w:sz w:val=\"{size_half_pt}\"/></w:rPr>"
        f"<w:t xml:space=\"preserve\">{escape(str(text))}</w:t></w:r>"
        "</w:p>"
    )


def _docx_table(rows: list[dict]) -> str:
    """Return an OpenXML table element from a list of row dicts."""
    if not rows:
        return _docx_paragraph("No table records available.")
    headers = list(rows[0].keys())
    parts = [
        "<w:tbl><w:tblPr><w:tblBorders>"
        "<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"0052CC\"/>"
        "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"0052CC\"/>"
        "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"E3E8F0\"/>"
        "</w:tblBorders></w:tblPr>"
    ]
    header_cells = "".join(
        f"<w:tc><w:p><w:r><w:rPr><w:b/><w:color w:val=\"0052CC\"/></w:rPr>"
        f"<w:t>{escape(str(h))}</w:t></w:r></w:p></w:tc>"
        for h in headers
    )
    parts.append(f"<w:tr>{header_cells}</w:tr>")
    for row in rows:
        cells = "".join(
            f"<w:tc><w:p><w:r><w:t>{escape(str(row.get(h, '')))}</w:t></w:r></w:p></w:tc>"
            for h in headers
        )
        parts.append(f"<w:tr>{cells}</w:tr>")
    parts.append("</w:tbl>")
    return "".join(parts)


def generate_word_report(
    project_title: str,
    parameters: dict,
    metrics: list[dict],
    table_rows: list[dict],
) -> bytes:
    """Compile project parameters, metrics, and table rows into a valid OpenXML Word (.docx) file."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    body_parts = [
        _docx_paragraph(f"{project_title} — Summary Report", bold=True, color="0052CC", size_half_pt=32),
        _docx_paragraph(f"Generated: {timestamp}", bold=False, color="5A6478", size_half_pt=18),
        _docx_paragraph("1. Project Parameters", bold=True, color="0A0F1C", size_half_pt=26),
    ]
    for key, val in parameters.items():
        body_parts.append(_docx_paragraph(f"• {key}: {val}"))

    body_parts.append(_docx_paragraph("2. Physicochemical & Summary Metrics", bold=True, color="0A0F1C", size_half_pt=26))
    for item in metrics:
        title = item.get("title", "Metric")
        value = item.get("value", "")
        rule = item.get("rule", "")
        status = item.get("status", "")
        body_parts.append(_docx_paragraph(f"• {title}: {value} ({rule} — {status})"))

    body_parts.append(_docx_paragraph("3. AutoDock Vina Calculated Scores Table", bold=True, color="0A0F1C", size_half_pt=26))
    body_parts.append(_docx_table(table_rows))

    document_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">"
        f"<w:body>{''.join(body_parts)}<w:sectPr/></w:body>"
        "</w:document>"
    )
    content_types_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/word/document.xml\" "
        "ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml\"/>"
        "</Types>"
    )
    rels_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" "
        "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" "
        "Target=\"word/document.xml\"/>"
        "</Relationships>"
    )

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as docx:
        docx.writestr("[Content_Types].xml", content_types_xml)
        docx.writestr("_rels/.rels", rels_xml)
        docx.writestr("word/document.xml", document_xml)
    buffer.seek(0)
    return buffer.getvalue()


def generate_combined_archive(
    project_title: str,
    parameters: dict,
    metrics: list[dict],
    table_rows: list[dict],
) -> bytes:
    """Bundle both the PDF and Word (.docx) reports into a single ZIP archive."""
    pdf_bytes = generate_pdf_report(project_title, parameters, metrics, table_rows)
    docx_bytes = generate_word_report(project_title, parameters, metrics, table_rows)

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("dragcorbs_report.pdf", pdf_bytes)
        zf.writestr("dragcorbs_report.docx", docx_bytes)
    buffer.seek(0)
    return buffer.getvalue()
