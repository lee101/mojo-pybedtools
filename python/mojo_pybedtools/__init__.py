"""Mojo-powered, in-memory interval arithmetic compatible with pybedtools' core API."""

from .bedtool import BedTool, Interval, create_interval_from_list, example_bedtool

__all__ = ["BedTool", "Interval", "create_interval_from_list", "example_bedtool"]
__version__ = "0.1.0"
