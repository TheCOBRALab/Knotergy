# Contribution Guide
This guide explains how to contribute to the project and keep our codebase clean, consistent, and maintainable.

## Debug

When developing code, **always build in Debug mode**. This enables runtime checks and detailed debug information, making it easier to catch potential issues early.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

## Style Guide

Consistent code style makes our codebase easier to read, maintain, and review. We use **clang-format** to automatically apply our project's style rules.

### Install clang-format

```bash
sudo apt install clang-format # Linux
brew install clang-format # macOS
```

### Format Your Code

To format all `.cpp` and `.hpp` files in the `src` & `test` directory:

```bash
cmake -S . -B build; cmake --build build --target format
```

> **Tip:** Run this before committing. Github CI will fail if you don't


## No Warnings Allowed

Your code **must compile without any warnings**. Treat every warning as an error.

### Why This Matters

Clean, warning-free builds:

* **Easier debugging:** Warning-free builds reduce noise when you're fixing issues.
* **Higher quality code:** Strict checks help avoid subtle bugs and undefined behavior.
* **Better design:** Many warnings are indicators of mismatched types or unsafe practices.
* **Keeps you safe:** If there's even one warning, ***I will find you.*** 

> **Tip:** If you find yourself repeatedly using `static_cast`, it may be a sign that you should reconsider your data types or design approach.

### Example

**❌ Problematic code:**

```c++
std::vector<int> numbers{-1, 2, 3};
for (int i = 0; i <= 2; ++i) {
    numbers[std::static_cast<std::size_t>(i)];
}
```

Here, the repeated `static_cast` is a red flag. You're forcing the compiler to silence a type mismatch instead of solving it at the root.

**✅ Better approach:**

```c++
std::vector<int> numbers{-1, 2, 3};
for (std::size_t i = 0; i < numbers.size(); ++i) {
    numbers[i];
}
```

By using the correct type (`std::size_t`) for indexing:
- The code is cleaner and easier to read.
- The variable matches the type expected by `std::vector::operator[]`.
- No signed/unsigned comparison warnings are generated.

> **Project Convention:** In Knotergy, `std::size_t` must be used for values representing indices. This makes intent clear and helps prevent bugs caused by using incorrect value types.


## How the program flows

main.cpp -> pipeline -> preprocessing -> loop_tree -> energy

main.cpp gets all the user inputs, and uses pipeline commands to parse those user inputs
those user inputs are then fed into the RNAProcessor in preprocessing which gathers data about the region (e.g. closed regions)
the data is then stored into processedRNAEntry.

We then feed the processedRNAEntry into the loop factory to create a tree of closed regions, where the root node is the entire sequence
and each node is a closed region. Nested closed regions are children of the parent closed region. LoopFactory also gathers data about each closed region.
These values are essential for the energy calculation.

The tree is then fed into the the energyCalculator named ComputeEnergy to give us our desired energy value