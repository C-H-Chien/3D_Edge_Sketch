# 3D Edge Sketch
### Research @ LEMS, Brown University

## Introduction:
This is an official repository for 3D Edge Sketch from multiple views with known poses. Basically, 3D Edge Sketch follows a few steps: _(i)_ detect [third-order edges](https://github.com/C-H-Chien/Third-Order-Edge-Detector) from each image, _(ii)_ pick two hypothesis views (with sufficient overlap in the image while maintaining some levels of baseline) and keep the rest of the views as validation views, _(iii)_ loop over all edges in the first hypothesis view, pair up with the edges falling in the epipolar wedge of the second hypothesis view, and seek support from the validation views, _(iv)_ triangulate all pairs of edges from the two hypothesis views to reconstruct 3D edges, _(v)_ reproject the 3D edges in order to pick the next two hypothesis views for the second round, _(vi)_ repeat steps _(iii)_ to _(vi)_ until almost all 2D edges claim support for 3D edges. <br />

The code is in a research/test status which renders most of the things unorganizned. In addition, steps _(v)_ and _(vi)_ are not included in the C++ code yet; they are manually done probably in the MATLAB code.

## Dependencies:
(1) CMake 3.14 or higher <br />
(2) Eigen 3.3.2 (higher version is not compatible) <br />

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
(3) For Brown University CCV Oscars users, manually add Eigen libraries cmake file to ``Eigen3_DIR``:
```bash
/gpfs/runtime/opt/eigen/3.3.2/share/eigen3/cmake/
```
keep pressing ``c`` for "[c] Configure" until "[g] Generate" appears. Press ``g`` to generate a make file. Finally, press ``q`` to exit the GUI.
(4) Compile the entire code:
```bash
make -j
```
and an executive file will be under ``/buid/bin``. <br />

## Running an Example
To get started, we use [ABC-NEF dataset](https://github.com/yunfan1202/NEF_code?tab=readme-ov-file#evergreen_treedataset) object 00000006. Third-order edges and absolute poses for all images are already processed which can be seen under ``datasets/ABC-NEF/00000006/Edges`` and ``datasets/ABC-NEF/00000006/RnT``, respectively. The macros in the ``definitions.h`` are set for this example, if otherwise specified. <br />
(1) If you are using multiple CPU cores, you can set the value of the macro ``NUM_OF_OPENMP_THREADS`` in ``definitions.h`` to the number of CPU cores you'd like to run 3D edge sketch in parallel. Remember to re-compile after the setting is done. <br />
(2) Simply run 3D Edge Sketch by
```bash
./edge_reconstruction-main
```
which is located under ``/build/bin``. <br />
(3) The output 3D edge points are written in the file under ``outputs/`` where the file is named by the name of the dataset, object, and the settings. <br />
(4) Visualize the 3D edge sketch using ``visualization/plot_3D_edge_sketch.m``. The 3D edges are under the coordinate of the first hypothesis view.

## Generating Third-Order Edges and Absolute Poses as Inputs for 3D Edge Sketch
- Third-order edges: use ``preprocesser/third_order_edge_detector/get_RO_Edges_List_in_dataset.m`` which helps generate ``Edge_*_t*.txt``. You may observe some ``.edg`` file which can be ignored for now as they are used by 3D curve sketch.
- Absolute Poses: use ``preprocesser/get_poses_from_ABC_NEF_dataset.m`` which reads dataset ground-truth file and transform all ground-truth poses to ``R_matrix.txt`` and ``T_matrix.txt``. The example code reads ``transforms_train.json`` file provided by the ABC-NEF dataset.

## TODOs
- [ ] Organize all messy MATLAB files 
- [ ] Organize the code to make it highly readable and run with no error
- [ ] Use YAML file to parse input arguments
- [ ] Implement multiple runs (automatically) to get a complete 3D edge sketch

## Contributors:
Yilin Zheng (yilin_zheng@alumni.brown.edu) <br />
Chiang-Heng Chien* (chiang-heng_chien@brown.edu) <br />
Qiwu Zhang (qiwu_zhang@brown.edu) <br />
*corresponding author


