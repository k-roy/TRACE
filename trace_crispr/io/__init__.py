"""
I/O modules for CRISPRo.

Author: Kevin R. Roy
"""

from .sample_key import (
    Sample,
    load_sample_key,
    create_sample_key_template,
)

from .output import (
    SampleResult,
    write_results_tsv,
    write_per_read_classifications,
    generate_summary_report,
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
