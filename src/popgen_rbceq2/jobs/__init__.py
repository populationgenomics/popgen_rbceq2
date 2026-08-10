"""Scripts a stage runs inside a batch job, invoked by path rather than imported.

A stage builds `python3 <this dir>/<name>_job.py --flag value ...` with
stage_support.build_python_command, so these modules run as scripts in the driver image. Keep
their imports to the parts of the package that carry no Hail Batch state.
"""
