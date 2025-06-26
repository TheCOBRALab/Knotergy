# Compute Energy Re-Write

> **Minimum‑free‑energy prediction of RNA secondary structures with pseudoknots**

Knotergy computes the free energy of a given structure



#### Steps for installation   
1. [Download the repository](https://github.com/TheCOBRALab/Knotergy) and extract the files onto your system.
2. From a command line in the root directory run
```
cmake -S . -B build
cmake --build build --parallel
``` 
If you want to run DEBUG mode, you can run
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```