# gssFortran

This repo is to test on using `class(*)` to bring arbitrary objective function into the algorithm and compare its efficiency with directly code golden section search for specific optimization.

My current conclusion is
- if you need to put `configurations` derived type into the `func_data` that the coding is doing, then the array in derived type should be in the heap memory rather than stack memory. In fortran syntax, the array should be declared using `allocatable` but not assign array `shape` when defined in the derived type.


