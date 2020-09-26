# Rat-FK3V-files
This repository consists of three folders:
* [FK3V](https://github.com/DGWhittaker/Rat-FK3V-files/tree/master/FK3V) contains C code implementing a modified Fenton-Karma model which reproduces rat ventricular APD, along with matplotlib script to compare APD restitution with experimental data:
<img src="https://github.com/DGWhittaker/Rat-FK3V-files/blob/master/FK3V/FK3V-fig.png">

* [Geometry-files](https://github.com/DGWhittaker/Rat-FK3V-files/tree/master/Geometry-files) contains 3D DTI-based geometry and myocyte orientation files in .txt format and a .vtk format geometry file; alpha and theta are rule-based implementations, whereas MH1, MH3, MH5 correspond to DTI1, DTI2, DTI3, respectively. Run `gcc eigenvec.c -lm -o eigenvec.o` followed by `./eigenvec.o` to compile and run code to extract eigenvector and myocyte orientation angle files for visualisation. An example figure showing myocyte orientation streamlines in a wedge is provided:
<img src="https://github.com/DGWhittaker/Rat-FK3V-files/blob/master/Geometry-files/streamlines-wedge.png" height="400">

* [Videos](https://github.com/DGWhittaker/Rat-FK3V-files/tree/master/Videos) contains supplementary videos from the paper

# Acknowledging this work

If you publish any work based on the contents of this repository please cite:

Whittaker, D. G., Benson, A. P., Teh, I., Schneider, J. E., Colman, M. A.
(2019).
[Investigation of the role of myocyte orientations in cardiac arrhythmia using image-based models](https://doi.org/10.1016/j.bpj.2019.09.041).
_Biophysical Journal_, 117, 2396-2408.