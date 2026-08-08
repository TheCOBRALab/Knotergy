# Knotergy: Pseudoknotted RNA Secondary Structure Energy Evaluator

> **Compute Gibbs free energy for RNA secondary structures, including pseudoknots**

[![GitHub release](https://img.shields.io/github/v/release/TheCOBRALab/Knotergy)](https://github.com/TheCOBRALab/Knotergy/releases)

Knotergy is a **C++ bioinformatics tool** for calculating the **Gibbs free energy (ΔG)** of RNA secondary structures, **including pseudoknots**.\
For non-pseudoknotted motifs, it uses the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) library.

## System Requirements

### Supported Operating Systems

- **Linux (64-bit)**
- **macOS (64-bit)**

## Install Knotergy via Conda

```
conda install cobralab::knotergy
conda activate
```

## Instructions to compile Knotergy locally

### Software Requirements

- **CMake** ≥ 3.15

---

## 1. Install CMake

Knotergy requires **CMake ≥ 3.15**.

If you have it installed, check your current version:

```bash
cmake --version
```

## How to install

### Linux

```bash
sudo apt update
sudo apt install cmake
```

### macOS

Install [Homebrew](https://brew.sh/) first, then:

```bash
brew install cmake
```

---

## 2. Download & Build Knotergy

### Clone the repository

```bash
git clone https://github.com/TheCOBRALab/Knotergy.git
cd Knotergy
```

Alternatively, download the zip:

```bash
curl -L -O https://github.com/TheCOBRALab/Knotergy/archive/refs/heads/main.zip
unzip main.zip
cd Knotergy-main
```

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release; cmake --build build --parallel
```

**Build Types:**

- **Release:** Most optimized
- **Debug:** Slower but catches hard to spot bugs (e.g. overflow)
- **Strict:** Same as Debug, but treats warnings as errors

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release; cmake --build build --parallel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug; cmake --build build --parallel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Strict; cmake --build build --parallel
```

**Build, Style, Test (one liner):**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug; cmake --build build --parallel; cmake --build build --target format; cd build; ctest --output-on-failure; cd ..
```

---

## 3. Usage

```bash
./build/Knotergy [flags]
```

Show all available options:

```bash
./build/Knotergy --help
```

---

## Examples

### Sequence + structure

```bash
./build/Knotergy -s AAAAUUU -r "((...))"
```

---

### Modified bases

```bash
# Use default modified base parameters
./build/Knotergy -s 6AAAUUU -r "((...))" -m

# Use a specific parameter file
./build/Knotergy -s 6AAAUUU -r "((...))" -m "./params/modified_bases/rna_mod_m6A_parameters.json"

# Use multiple parameter files
./build/Knotergy -s 6AAAUUU -r "((...))" \
  -m "./params/modified_bases/rna_mod_m6A_parameters.json" \
  -m "./params/modified_bases/rna_mod_pseudouridine_parameters.json"

# Use all parameter files in a directory
./build/Knotergy -s 6AAAUUU -r "((...))" -m "./params/modified_bases/"
```

---

### Input / output

```bash
# Read from input file
./build/Knotergy -i ./templates/sample_input.fa

# Write output to file
./build/Knotergy -i ./templates/sample_input.fa > output.txt
```

---

### Pseudoknot energy rounding

Use `-e` or `--round` to control rounding of **individual pseudoknot stacking and internal-loop energy contributions after their pseudoknot multipliers are applied**. Rounding occurs before these contributions are added to the total energy. Other energy terms are not affected.

```bash
./build/Knotergy -i input.fa -e        # Banker's rounding
./build/Knotergy -i input.fa -e 2      # Round to nearest
./build/Knotergy -i input.fa -e2       # Compact form
./build/Knotergy -i input.fa --round 3
```

| Value | Method                                         |
| ----: | ---------------------------------------------- |
|   `0` | No rounding *(used by HotKnots)*               |
|   `1` | Banker's rounding *(used by HFold and Knotty)* |
|   `2` | Round to nearest integer                       |
|   `3` | Round down (`floor`)                           |
|   `4` | Round up (`ceil`)                              |
|   `5` | Truncate toward zero                           |


Without `-e`, no rounding is applied. Using `-e` without a value selects Banker's rounding.

---

## Command-line options

### 🔧 General options

- `-h`, `--help`  
  Show this help message and exit

- `-V`, `--version`  
  Print version information and exit

- `-v`, `--verbose`  
  Enable detailed output, including per-motif energy breakdown

---

### 📥 Input options

- `-s`, `--sequence`  
  RNA sequence input

- `-r`, `--structure`  
  RNA secondary structure (dot-bracket notation)

- `-i`, `--input`  
  Input file in FASTA format (see `templates/sample_input.fa`)

---

### ⚙️ Parameter files

- `-P`, `--paramFile`  
  ViennaRNA parameter file

- `-k`, `--pk-paramFile`  
  Pseudoknot parameter file

- `-m`, `--mod-params`  
  Modified base parameters (ViennaRNA format; accepts file(s) or directory)

---

### 🧪 Calculation options

* `-e`, `--round [0-5]`
  Select how pseudoknot stacking and internal-loop energy contributions are rounded after applying their pseudoknot multipliers. Use `0` for no rounding, `1` for Banker's rounding, or `2`–`5` for the other supported methods. (default 0)


- `-d`, `--dangle [0-3]`  
  Dangle model to use (default: `2`)
