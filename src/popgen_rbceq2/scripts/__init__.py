"""Developer scripts, run by hand rather than by the workflow.

gen_bg_resources.py regenerates the committed `resources/bg_*.<genome>.*` from the rbceq2
allele database; bg_db.py is the parsing it shares with the QC job. gen_off_design_sites.py
subtracts one exome capture design from those defining sites and commits the result beside
them, with a manifest row saying what it was built from.
"""
