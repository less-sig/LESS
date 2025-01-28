Memory Analysis:
================

Run `python bench.py` to start the stack analysis. Standardly the AVX2 
implementation will be analyzed. To change you need to change line 20-22
in `CMakeLists.txt` to chose a different implementation.

The bash script `bench.sh` is automatically building the C library and running
`valgring` to get the performance data. Note: the bash script cannot analyse 
those data. This can only be done by the python script.
