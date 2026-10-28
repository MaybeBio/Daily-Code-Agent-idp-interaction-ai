#!/usr/bin/env python3
"""
generate_report.py
Combines PNG plots into a single PDF report.

Usage:
    python generate_report.py <dataset_name> <output_pdf> <plot1.png> [plot2.png ...]
"""

import sys
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, HRFlowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from PIL import Image as PILImage


# Map filename fragments to human-readable section headings (order matters —
# more specific fragments must come before more general ones)
PLOT_LABELS = [
    ("plddt_density_scipy",                          "pLDDT Distribution (SciPy)"),
    ("plddt_gmm_scatter_disorder",                   "Proteome Structural Confidence Landscape — coloured by mean disorder"),
    ("plddt_gmm_scatter",                            "Proteome Structural Confidence Landscape"),
    ("metapredict_mean_disorder_density_scipy",      "Metapredict — Mean Disorder Density (SciPy)"),
    ("metapredict_mean_disorder_density_statsmodels","Metapredict — Mean Disorder Density (Statsmodels)"),
    ("metapredict_mean_disorder_hist",               "Metapredict — Mean Disorder Histogram"),
    ("metapredict",                                  "Metapredict Disorder"),
    ("psauron",                                      "pSAURON Distribution"),
    ("interpro_merged_coverage_distribution",        "InterPro Domain Coverage — Per-Protein Distribution"),
    ("interpro_domain_coverage",                     "InterPro Domain Coverage — Per-Library"),
]


def label_for(png_path):
    name = os.path.basename(png_path).lower()
    for fragment, label in PLOT_LABELS:
        if fragment in name:
            return label
    return os.path.splitext(os.path.basename(png_path))[0]


def fit_image(path, max_width, max_height):
    with PILImage.open(path) as im:
        native_w, native_h = im.size
    scale = min(max_width / native_w, max_height / native_h)
    return Image(path, width=native_w * scale, height=native_h * scale)


def build_report(dataset_name, output_pdf, png_files):
    page_width, page_height = A4
    margin = 2 * cm
    img_width = page_width - 2 * margin
    img_max_height = page_height * 0.72

    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=A4,
        leftMargin=margin,
        rightMargin=margin,
        topMargin=margin,
        bottomMargin=margin,
        title=f"{dataset_name}",
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "ReportTitle",
        parent=styles["Title"],
        fontSize=17,
        leading=22,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#1a1a2e"),
        spaceAfter=6,
    )

    section_style = ParagraphStyle(
        "SectionHeading",
        parent=styles["Heading2"],
        fontSize=12,
        leading=16,
        textColor=colors.HexColor("#16213e"),
        spaceBefore=10,
        spaceAfter=6,
    )

    story = []

    # Title + decorative rule
    story.append(Spacer(1, 0.6 * cm))
    story.append(Paragraph(
        f"{dataset_name}",
        title_style
    ))
    story.append(Spacer(1, 0.3 * cm))
    story.append(HRFlowable(
        width="100%", thickness=1.5,
        color=colors.HexColor("#4a4e69"), spaceAfter=0.6 * cm
    ))

    for i, png in enumerate(png_files):
        story.append(Paragraph(label_for(png), section_style))
        story.append(fit_image(png, img_width, img_max_height))
        if i < len(png_files) - 1:
            story.append(PageBreak())

    doc.build(story)
    print(f"Report written to: {output_pdf}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: generate_report.py <dataset_name> <output_pdf> <plot1.png> [plot2.png ...]"
        )
    dataset_name = sys.argv[1]
    output_pdf   = sys.argv[2]
    png_files    = sys.argv[3:]
    build_report(dataset_name, output_pdf, png_files)


