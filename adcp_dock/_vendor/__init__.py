"""Modules shared with the HOPE peptide pipeline.

This tool prepares a receptor, builds a grid target, docks with ADCP and
rescores with MM-GBSA using the same code as that pipeline's hope_core
package, copied here without code changes so that a fix in one reaches the
other. The copy is kept apart under _vendor/ so the boundary is visible.

Not copied: pose pairing and site clustering (pairs.py, sites.py), which
belong to a growth loop this tool does not run, and tools.py, since
adcp_dock/tools.py also finds the Amber programs. estimate.py is an
excerpt.

scripts/refresh_vendor.py refreshes the copy and reports code that differs.
"""
