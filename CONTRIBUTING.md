# Contribution Guide

## Style Guide
How to auto style your code using out repo's style guide

``` bash
sudo apt install clang-format # Linux
brew install clang-format     # MacOS

find ./src \( -iname "*.cpp" -o -iname "*.hpp" \) -print0 | xargs -0 clang-format -i # Formats your code
```

## How the program flows
main.cpp -> pipeline -> preprocessing -> loop_tree -> energy

main.cpp gets all the user inputs, and uses pipeline commands to parse those user inputs
those user inputs are then fed into the RNAProcessor in preprocessing which gathers data about the region (e.g. closed regions)
the data is then stored into processedRNAEntry.

We then feed the processedRNAEntry into the loop factory to create a tree of closed regions, where the root node is the entire sequence
and each node is a closed region. Nested closed regions are children of the parent closed region. LoopFactory also gathers data about each closed region.
These values are essential for the energy calculation.

The tree is then fed into the the energyCalculator named ComputeEnergy to give us our desired energy value