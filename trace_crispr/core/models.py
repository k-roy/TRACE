"""
Data models for TRACE CRISPR analysis.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class EditingTemplate:
    """
    Represents a template for CRISPR-mediated editing.

    Attributes:
        template_id: Unique identifier for the template
        sequence: DNA sequence of the template
        expected_edit: Expected edit outcome (e.g., "C>T at position 123")
        guide_sequence: Guide RNA sequence used
        description: Human-readable description
    """
    template_id: str
    sequence: str
    expected_edit: Optional[str] = None
    guide_sequence: Optional[str] = None
    description: Optional[str] = None

    def __repr__(self) -> str:
        return f"EditingTemplate(id={self.template_id}, edit={self.expected_edit})"
