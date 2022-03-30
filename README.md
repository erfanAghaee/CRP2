CRP2.0
======================================
CRP2.0 (Co-operation between Routing and Placement) is a VLSI global routing with cell movement tool developed Erfan Aghaeekiasaraeee in the research team supervised by Prof. Laleh Behjat at University of Calgary. This source code is evolved on CUGR (https://github.com/cuhk-eda/cu-gr)

## 1. How to Build

**Step 1:** Download the source code. For example,
```bash
$ git clone https://github.com/erfanAghaee/CRP2.git
```

**Step 2:** Go to the project root and build by
```bash
$ cd CRP2
$ ./build.sh #scripts/build.py -o release
$ cd ../buildCRP2/
$ ./make.sh 
$ cd ../buildTriton/
$ cmake ../CRP2/triton/
$ ./make.sh 
$ cd run/
$ ./run_exps

Note that this will generate two folders under the root, `build` and `run` (`build` contains intermediate files for build/compilation, while `run` contains binaries and auxiliary files).
More details are in [`scripts/build.py`](scripts/build.py).

### 1.1. Dependencies

* [GCC](https://gcc.gnu.org/) (version >= 5.5.0) or other working c++ compliers
* [CMake](https://cmake.org/) (version >= 2.8)
* [Boost](https://www.boost.org/) (version >= 1.58)
* [Python](https://www.python.org/) (version 3, optional, for utility scripts)
* [Innovus®](https://www.cadence.com/content/cadence-www/global/en_US/home/tools/digital-design-and-signoff/soc-implementation-and-floorplanning/innovus-implementation-system.html) (version 18.1, optional, for design rule checking and evaluation)
* [Rsyn](https://github.com/RsynTeam/rsyn-x) (a trimmed version is used, already added under folder `rsyn`)

