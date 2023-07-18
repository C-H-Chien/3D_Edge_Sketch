# Edge-Based Reconstruction
### Research @ LEMS, Brown University (CVPR 2022)

## Dependencies:
(1) CMake 3.14 or higher <br />
(2) Eigen 3.3.2 or higher <br />

## How to build and compile the code
(1) After cloning the repo, cd to the repo folder and create a 'build' directory and enter it
```bash
mkdir build && cd build
```
(2) Create a make file <br />
```bash
ccmake ..
```
ans you will see a GUI where you can explicitly provide directories. <br />
(3) For CCV users, manually add Eigen libraries cmake file to ``Eigen3_DIR`` then press ``q`` to exit the GUI:
```bash
/gpfs/runtime/opt/eigen/3.3.2/share/eigen3/cmake/
```
(4) Compile the entire code:
```bash
make -j
```
(5) An executive file will be under ``/buid/bin``.
(6) CCV users, remember to launch an interactive CPU node to run the code:
```bash
interact -n 8 -t 03:00:00 -m 3g
```
which means, an interactive node with 8-cores CPU, 3 hours of usage, and holding 3G memories.