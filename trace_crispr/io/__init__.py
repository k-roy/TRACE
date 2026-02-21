"""
I/O modules for CRISPRo.

Author: Kevin R. Roy
"""

from .output import (
    SampleResult,
    generate_summary_report,
    write_per_read_classifications,
    write_results_tsv,
)
from .sample_key import (
    Sample,
    create_sample_key_template,
    load_sample_key,
)

__all__ = [
    'Sample',
    'load_sample_key',
    'create_sample_key_template',
    'SampleResult',
    'write_results_tsv',
    'write_per_read_classifications',
    'generate_summary_report',
]
