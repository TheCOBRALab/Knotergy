# Compute Energy Re-Write

> **Minimum‑free‑energy prediction of RNA secondary structures with pseudoknots**

Knotergy computes the free energy of a given structure
## System Requirements
HFold needs a 64bit Linux or macOS operating system to run

### Software Requirements:
- CMake 3.15 or greater
- ViennaRNA 2.7.0

### ViennaRNA package installation
You will need to install ViennaRNA package in order to compile HFold

1. Download the ViennaRNA package from [ViennaRNA-2.7.0.tar.gz](https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz)
or using the command-line 
```bash
curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz
```


2. Install the ViennaRNA package:
```bash
tar -zxvf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0
./configure --without-perl
sudo make -j$(nproc)              # Linux
sudo make -j$(sysctl -n hw.ncpu)  # macOS
sudo make install
```

#### Common issues:
1. No admin access: Build it into a directory where you do have access. (prefix can be any directory you want)
```bash
mkdir -p ~/local
./configure --without-perl --prefix=$HOME/local
make -j$(nproc)              # Linux
make -j$(sysctl -n hw.ncpu)  # macOS
make install
```

2. Template error in dlib: Fix the bug in the code using your command-line (usually happens on newer macs)
```bash
sed -i.bak 's/::template go(/::template go<>(/' src/dlib-19.24/dlib/global_optimization/find_max_global.h
```

For more details about ViennaRNA, see https://github.com/ViennaRNA/ViennaRNA  

### CMake installation


[CMake](https://cmake.org/install/) version 3.15 or greater must be installed in a way that HFold can find it.    

To check if your system already has CMake, run this in terminal:     
```
cmake --version
```

#### Linux:
```
sudo apt update
sudo apt install cmake
```

#### Mac:    
Homebrew is required to download CMake on Mac.
Run this command to install homebrew. 
```  
curl -L -O https://github.com/TheCOBRALab/Knotergy/archive/refs/heads/main.zip
``` 
You may be prompted to add homebrew to PATH. If prompted, please follow the instructions presented on the terminal.

Once finished, install CMake
```   
brew install cmake   
``` 

## Install Knotergy
1. [Download the repository](https://github.com/TheCOBRALab/HFold) either manually or using the command-line 
``` bash
# using git
git clone https://github.com/TheCOBRALab/Knotergy.git

# using curl
/bin/bash -c "$(curl -fsSL https://github.com/TheCOBRALab/Knotergy/archive/refs/heads/main.zip)"
```

2. Unzip the file and enter the root dir of knotergy (where README.md is)
``` bash
unzip Knotergy-main.zip
cd Knotergy-main
```

3. Build the program
```bash
cmake -S . -B build
cmake --build build --parallel
```

If you want to run DEBUG mode, you can run

``` bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```


## Usage

```
Usage: ./build/Knotergy [flags]
```
