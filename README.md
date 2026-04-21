# Knotergy: Pseudoknotted RNA Secondary Structure Energy Evaluator

> **Compute Gibbs free energy for RNA secondary structures, including pseudoknots**

[![GitHub release](https://img.shields.io/github/v/release/TheCOBRALab/Knotergy)](https://github.com/TheCOBRALab/Knotergy/releases)

Knotergy is a **C++ bioinformatics tool** for calculating the **Gibbs free energy (ΔG)** of RNA secondary structures, **including pseudoknots**.\
For non-pseudoknotted motifs, it uses the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) library.

## System Requirements

### Supported Operating Systems

- **Linux (64-bit)**
- **macOS (64-bit)**

### Software Requirements

- **CMake** ≥ 3.15
- **ViennaRNA** ≥ 2.7.0

---

## 1. Install ViennaRNA

Knotergy requires ViennaRNA to compile and run.

### Option 1: Install via Conda (Easiest)

1. **Install ViennaRNA through conda, and activate your conda environment**

   ```bash
   conda install -c bioconda viennarna
   conda activate
   ```

### Option 2: Manual install (Advanced / for development)

1. **Download ViennaRNA 2.7.2:**

   ```bash
   curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.2/ViennaRNA-2.7.2.tar.gz
   ```

2. **Install:**

   ```bash
   tar -zxvf ViennaRNA-2.7.2.tar.gz
   cd ViennaRNA-2.7.2
   ./configure --without-perl
   make -j$(nproc)              # Linux
   make -j$(sysctl -n hw.ncpu)  # macOS
   sudo make install
   ```

### Common Installation Issues

- **No admin access:**

  ```bash
  mkdir -p ~/local
  ./configure --without-perl --prefix=$HOME/local
  make -j$(nproc)              # Linux
  make -j$(sysctl -n hw.ncpu)  # macOS
  make install
  ```

- **macOS: “not 8-byte aligned”**

  ```bash
  make clean
  export AR=/usr/bin/ar
  export RANLIB=/usr/bin/ranlib
  ./configure --without-perl
  make -j$(sysctl -n hw.ncpu)
  sudo make install
  ```

For full details, see the [ViennaRNA GitHub repo](https://github.com/ViennaRNA/ViennaRNA).

---

## 2. Install CMake

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

## 3. Download & Build Knotergy

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

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release; cmake --build build --parallel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug;   cmake --build build --parallel
```

**Build & Test (one liner):**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug; cmake --build build --parallel; cd build; ctest --output-on-failure; cd ..
```

---

## 4. Usage

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

- `-p`, `--paramFile`  
  ViennaRNA parameter file

- `-k`, `--pk-paramFile`  
  Pseudoknot parameter file

- `-m`, `--mod-params`  
  Modified base parameters (ViennaRNA format; accepts file(s) or directory)

---

### 🧪 Calculation options

- `-e`, `--round`  
  Round pseudoknot energies to integers. This matches the behavior of some predictors such as HFold and Knotty, which trade accuracy for speed. Knotergy itself does not gain performance from this option.

- `-d`, `--dangle`  
  Dangle model to use (default: `2`)
