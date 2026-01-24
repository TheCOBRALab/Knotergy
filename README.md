# Knotergy: Pseudoknotted RNA Secondary Structure Energy Evaluator

> **Compute Gibbs free energy for RNA secondary structures, including pseudoknots**

[![GitHub release](https://img.shields.io/github/v/release/TheCOBRALab/Knotergy)](https://github.com/TheCOBRALab/Knotergy/releases)


Knotergy is a **C++ bioinformatics tool** for calculating the **Gibbs free energy (ΔG)** of RNA secondary structures, **including pseudoknots**.\
It uses the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) library for non-pseudoknotted calculations.
## System Requirements

### Supported Operating Systems

* **Linux (64-bit)**
* **macOS (64-bit)**

### Software Requirements

* **CMake** ≥ 3.15
* **ViennaRNA** 2.7.0 & 2.7.1
---

## 1. Install ViennaRNA

Knotergy requires ViennaRNA to compile and run.

1. **Download ViennaRNA 2.7.1:**

   ```bash
   curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.1/ViennaRNA-2.7.1.tar.gz
   ```

2. **Install:**

   ```bash
   tar -zxvf ViennaRNA-2.7.1.tar.gz
   cd ViennaRNA-2.7.1
   ./configure --without-perl
   sudo make -j$(nproc)              # Linux
   sudo make -j$(sysctl -n hw.ncpu)  # macOS
   sudo make install
   ```

### Common Installation Issues

* **No admin access:**

  ```bash
  mkdir -p ~/local
  ./configure --without-perl --prefix=$HOME/local
  make -j$(nproc)              # Linux
  make -j$(sysctl -n hw.ncpu)  # macOS
  make install
  ```

For full details, see the [ViennaRNA GitHub repo](https://github.com/ViennaRNA/ViennaRNA).

---

## 2. Install CMake

Knotergy requires **CMake ≥ 3.15**.

Check your current version:

```bash
cmake --version
```

### Linux:

```bash
sudo apt update
sudo apt install cmake
```

### macOS:

Install [Homebrew](https://brew.sh/) first, then:

```bash
brew install cmake
```

---

## 3. Install Knotergy

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
cmake -S . -B build
cmake --build build --parallel
```

**Debug build:**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
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

Run `--help` for available options:

```bash
./build/Knotergy --help
```