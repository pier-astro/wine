"""Read and evaluate standard XSPEC table models without XSPEC."""

from .table import TableSpectrum, XspecTableModel, read_table_model

__all__ = ["TableSpectrum", "XspecTableModel", "read_table_model"]
