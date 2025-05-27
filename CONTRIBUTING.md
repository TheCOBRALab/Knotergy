# Contribution Guide

## Style Guide
How to auto style your code using out repo's style guide

``` bash
sudo apt install clang-format # Linux
brew install clang-format     # MacOS

find ./src \( -iname "*.cpp" -o -iname "*.hpp" \) | xargs clang-format -i # Formats your code
```