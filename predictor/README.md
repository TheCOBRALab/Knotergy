# Proof of concept predictor

Very simple predictor that uses Knotergy to help with its predictions. 
Can be tweaked in `src/settings.py`

DISCLAIMER: This predictor's results are not reliable and can produce unrealistic structures with favorable energies. This is purely meant as a proof of concept and is NOT to be used in real world applications.

## How to use

### Install
Either install Knotergy locally using conda
```
conda install cobralab::knotergy
conda activate
```

or build the code manually

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release; cmake --build build --parallel
```

### Run
Assuming you're in the `/Predictor` directory, run the following command

If installed via conda:
```
python3 main.py -s "seq"
```

> **Note:** `-i`, `-s`, `-e`, `-P`, and other flags SHOULD work. (Emphasis on should. I barely tested this)


If compiled locally:
```
python3 main.py --knotergy ../build/Knotergy -s "seq"
```

> **Note:** The `--knotergy flag` is used to set the path to the executable
