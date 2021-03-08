# WoodPixel

Supplementary code for Computational Parquetry: Fabricated Style Transfer with Wood Pixels, ACM Transactions on Graphics 39(2), 2020

### Authors
* Julian Iseringhausen, <[opensource@iseringhausen.graphics](mailto:opensource@iseringhausen.graphics)>
* Matthias Hullin, University of Bonn, <[hullin@cs.uni-bonn.de](mailto:hullin@cs.uni-bonn.de)>

### Build instructions (Ubuntu 20.04.1)
The following instructions are tested on a minimal Ubuntu 20.04.1 installation with gcc 9.3.0.

Install necessary packages.
```
sudo apt update && sudo apt install \
build-essential \
cmake \
git \
libboost-all-dev \
libopencv-dev \
```

Configure and build CMake project.
```
mkdir build
cd build
cmake ..
make
```

### Links
[Paper](https://light.cs.uni-bonn.de/computational-parquetry-fabricated-style-transfer-with-wood-pixels/)
