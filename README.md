# Multiview Edge-Based Reconstruction
### Research @ LEMS, Brown University

This is a branch where GPU implementation of multiview edge-based reconstruction is, which will be updated by following the updates from the main branch.

## Dependencies:
(1) CMake 3.14 or higher <br />
(2) gcc 10.2 or higher (lower might also doable but not tested yet) <br />
(3) Eigen 3.3.2 or higher <br />
(4) CUDA 11.1.1 or higher (if NVIDIA Nsight profiler is used, cuda/12.2.2 must be used)

## How to build and compile the code
(1) After cloning the repo, cd to the repo folder and create a 'build' directory and enter it
```bash
mkdir build && cd build
```
(2) Create a make file <br />
```bash
ccmake ..
```
and you will see a GUI where you can explicitly provide directories. <br />
(3) For Brown University CCV Oscars users, manually add Eigen libraries cmake file to ``Eigen3_DIR`` then press ``q`` to exit the GUI:
```bash
/gpfs/runtime/opt/eigen/3.3.2/share/eigen3/cmake/
```
(4) Compile the entire code:
```bash
make -j
```
(5) An executive file will be under ``/buid/bin``. <br />
(6) Brown University CCV Oscars users, remember to launch an interactive GPU node before running the code:
```bash
interact -q gpu -t 03:00:00 -g 1 -n 4
```
which means, an interactive GPU node with 1 GPU device, 4-cores CPU, and 3 hours of usage.

## How to profile using NVIDIA Nsight
Simply use
```bash
nsys profile ./edge_reconstruction-main
```
and you will get a ``report1.nsys-rep`` report which can be opened and visualized by the Nsight profiler tool.


