"""
aKmerBroom: Ancient oral DNA decontamination using Bloom filters on k-mer sets

This package provides tools for removing ancient oral DNA contamination from ancient DNA samples
using a two-step algorithm based on k-mer analysis.
"""

__version__ = "1.0.0"
__author__ = "Camila Duitama González"

from .main import main

__all__ = ["main"]