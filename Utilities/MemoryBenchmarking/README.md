Memory Analysis:
================

Run `python bench.py` to start the stack analysis. Standardly the AVX2 
implementation will be analyzed. To change you need to change line 20-22
in `CMakeLists.txt` to chose a different implementation.

The bash script `bench.sh` is automatically building the C library and running
`valgring` to get the max stack size.


Result:
======

Benchmarked on a Ryzen 7600X:

``` python
# reference impl
{1: {'INTERMEDIATE': {'keygen': 82472, 'sign': 171736, 'verify': 110216}, 'SHORT_SIG': {'keygen': 82616, 'sign': 159992, 'verify': 110104}, 'BALANCED': {'keygen': 82376, 'sign': 172440, 'verify': 110232}}, 3: {'INTERMEDIATE': {'keygen': 0, 'sign': 0, 'verify': 0}, 'SHORT_SIG': {'keygen': 204840, 'sign': 836936, 'verify': 300120}, 'BALANCED': {'keygen': 204760, 'sign': 754728, 'verify': 299864}}, 5: {'INTERMEDIATE': {'keygen': 0, 'sign': 0, 'verify': 0}, 'SHORT_SIG': {'keygen': 381496, 'sign': 1270320, 'verify': 525104}, 'BALANCED': {'keygen': 381400, 'sign': 1702640, 'verify': 591816}}}

# AVX2 impl.
{1: {'INTERMEDIATE': {'keygen': 119496, 'sign': 206328, 'verify': 142152}, 'SHORT_SIG': {'keygen': 119640, 'sign': 194584, 'verify': 142008}, 'BALANCED': {'keygen': 119384, 'sign': 207032, 'verify': 142136}}, 3: {'INTERMEDIATE': {'keygen': 0, 'sign': 0, 'verify': 0}, 'SHORT_SIG': {'keygen': 284872, 'sign': 977752, 'verify': 375480}, 'BALANCED': {'keygen': 284808, 'sign': 886072, 'verify': 375224}}, 5: {'INTERMEDIATE': {'keygen': 0, 'sign': 0, 'verify': 0}, 'SHORT_SIG': {'keygen': 477128, 'sign': 1403800, 'verify': 618520}, 'BALANCED': {'keygen': 477032, 'sign': 1854648, 'verify': 684952}}}
```
