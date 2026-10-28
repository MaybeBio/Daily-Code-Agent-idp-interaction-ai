"""Reproducible supplementary-figure builders."""

from pathlib import Path

from .shared import DEFAULT_OUT, setup_style, write_readme
from .figure_s01 import figure_s1_density
from .figure_s02 import figure_s2_contact_maps
from .figure_s03 import figure_s3_compensation
from .figure_s06 import figure_s6_tica
from .figure_s07 import figure_s7_its
from .figure_s08 import figure_s8_ck


FIGURE_BUILDERS = {
    1: figure_s1_density,
    3: figure_s2_contact_maps,
    4: figure_s6_tica,
    5: figure_s7_its,
    6: figure_s8_ck,
    8: figure_s3_compensation,
}


def run_figure(number: int, output_dir: Path = DEFAULT_OUT) -> None:
    """Generate one supplementary figure and refresh the output README."""
    if number not in FIGURE_BUILDERS:
        raise ValueError(f"Unknown supplementary figure: S{number}")
    setup_style()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    FIGURE_BUILDERS[number](output_dir)
    write_readme(output_dir)
