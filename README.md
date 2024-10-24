# Open Source: Geophysical 3D/2D Finite Difference Modeling and Marchenko Algorithms

This repository contains tools for geophysical modeling, including:
- 3D/2D Finite Difference (FD) modeling
- Marchenko algorithms
- 2D/3D x-w migration utilities

These tools are essential for seismic exploration, imaging, and multiple elimination tasks in inhomogeneous media. The code supports seismic interferometry and other advanced geophysical techniques.

Table of Contents
---------
- [Acknowledgments](#acknowledgement)
- [Features](#features)
- [Installation](#installation)
- [Reproducing](#reproducing)
- Demos:
    - [Finite Difference Modeling: FDELMODC](#finite-difference-modeling-fdelmodc)
    - [Marchenko Method: Marchenko](#marchenko-method--marchenko)
- [References](#references)
- [License](#license)

- Additional Installations
    - [Multidimensional Deconvolution](#mdd)
    - [Math Kernel Library](#mkl-libraries)
    - [Seismic Unix](#su)
    - [zfp](#zfp)
    - [FDACRTMC](#fdacrtmc)
- [Updates and Latest Version](#updates-and-latest-version)

ACKNOWLEDGEMENT
---------
This work received funding from the European Research Council (grant 742703) and the  NWO Domain Applied and Engineering Sciences (grant 13939).

FEATURES
---------
- **Finite Difference (FD) Modeling:** Includes 2D and 3D FD models for seismic interferometry, simulation, and analysis.
- **Marchenko Method:** Implements advanced algorithms for seismic wavefield extrapolation and multiple elimination.
- **x-w Migration Utilities:** Provides tools for x-w (space-frequency) migration in 2D/3D.
- **Seismic Unix Compatibility:** Integration with routines from Seismic Unix.

INSTALLATION
-------------
See the separate INSTALL file.


REPRODUCING
----------
Almost all Figures in the papers mentioned above can be reproduced with the software in this repository. Please see the file REPRODUCE for further instructions



Finite Difference Modeling: FDELMODC
------------------------------------
To run a demo for finite difference modeling, navigate to the fdelmodc/demo directory and execute the following script:



The demo directory contains many scripts which demonstrate the different possibilities of the modeling program.

To reproduce the Figures shown in the GEOPHYSICS manuscript "Finite-difference modeling experiments for seismic interferometry" the scripts in FiguresPaper directory can be used. Please read the README in the FiguresPaper directory for more instructions and guidelines.

An extensive manual of fdelmodc can be found in doc/fdelmodcManual.pdf

Marchenko method : MARCHENKO
----------------------------

To run the demos for the Marchenko method, navigate to the bin/marchenko directory and run the scripts found in:

marchenko/demo/oneD

marchenko/demo/twoD


To reproduce the Figures shown in the GEOPHYSICS paper "Implementation of the Marchenko method" the scripts in marchenko/demo/oneD directory can be used. The README in this directory gives more instructions and guidelines.

To reproduce the Figures shown in the Scientific Reports paper "Virtual acoustics in inhomogeneous media with single-sided access" the scripts in marchenko/demo/ScientificReports directory can be used. The README in this directory gives more instructions and guidelines.

To reproduce the Figures shown in the GEOPHYSICS paper "Implementation of the Marchenko Multiple Elimination algorithm" the scripts in marchenko/demo/mme directory can be used. The README_PRIMARIES in this directory gives more instructions and guidelines.

A brief manual about the MME program 'marchenko_primaries' can be found in doc/MMEmanual.pdf

A more extensive manual about Marchenko can be found in doc/marchenkoManual.pdf


MDD
---
Multidimensional deconvolution is a geophysical technique used for multidimensional wavefield reconstruction. The MDD kernels depend on **BLAS** and **LAPACK** calls. Free downloads of these libraries can be found on


https://www.netlib.org/blas/index.html
https://www.netlib.org/lapack/index.html


MDD with **LSQR** uses the Fortran90 zLSQR code (for real or complex A) from the Systems Optimization Laboratory of Stanford University


https://stanford.edu/group/SOL/software/lsqr/

MKL libraries
-------------

The Math Kernel Library (MKL) is a set of mathematical routines designed to accelerate performance on Intel processors. It helps optimize performance for x84_64 architectures.

If you are running on x64_86 processors you can download (for free) Intel's highly optimised MKL package:


https://software.intel.com/en-us/mkl/choose-download
 
These Libraries include highly optimised libraries of **BLAS**, **LAPACK**, and **FFT(W)**.
Usually the MKL libraries are installed in $MKLROOT. If that variable is not set, you can try to find the correct path by searching for one of the libraries:


find /opt/intel -name libmkl_gf_lp64.so


and adjust MKLROOT in Make_include accordingly.
You can also completely disable the use of MKL by commenting out the MKL parts in Make_include.


In case MKL is installed on your system there is no need to install the netlib packages mentioned for MDD.


SEISMIC UNIX
-----------

The output file format is based on Seismic Unix (SU) with a file extension '.su'. These files can be opened and processed with Seismic Unix Software.

If you want to use the .su files with SU from CWP:
git clone https://github.com/JohnWStockwellJr/SeisUnix


==> Please make sure that SU is compiled without XDR (in $CWPROOT/Makefile.config make sure that XDRFLAG is NOT enabled). The SU output files of fdelmodc are all based on local IEEE data.
To exclude the XDRFLAG in SU you have to use the following line in $CWPROOT/src/Makefile.config around line 35:


XDRFLAG =


so the XDRFLAG is empty. If you have made that change you have to remake SU with the commands:
make remake
make xtremake




ZFP
---
The fdacrtmc and the 3D Marchenko code makes use of ZFP compression to store snapshots in CPU memory or do efficient IO. This package is included in this repository for your convenience. The latest package and detailed explanation can be found on:



https://github.com/LLNL/zfp


and written by Peter Lindstrom.



FDACRTMC
--------
fdacrtmc uses **FFTW** and the wisdom computations are stored on disk for re-usage.  This directory is defined in fdacrtmc.h


#ifndef WISDOMDIR
#define WISDOMDIR "/tmp/fftw/"
#endif


For the code to run properly the directory /tmp/fftw/ must exist. It is also possible to compile the fdacrtmc code that defines
WISDOMDIR. In fdacrtmc/Makefile add:


CFLAGS  += -DWISDOMDIR="/directory/that/exists"


or you can change the name of WISDOMDIR in fdacrtmc.h


The Finite Different based RTM code (FDACRTMC) will not be installed if MKL is not available. This code depends on FFTW and we do want to include the full FFTW source code package to this. Our aim is to keep the compilation of the code in the GitLab repository as simple as possible with the smallest number of external libraries.


LICENSE
---------
THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THIS COMMON PUBLIC LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT.


A copy of this license can be found in the file 'Common_Public_License.txt' in the directory where you have found this README.


http://www.opensource.org/licenses/cpl1.0.php


SU
--
Some routines are from Seismic Unix and include the SU LEGAL_STATEMENT in the source code.


Copyright (c) 2017 by the Society of Exploration Geophysicists.
For more information, go to http://software.seg.org/2017/00XX .
You must read and accept usage terms at:
http://software.seg.org/disclaimer.txt before use.


REFERENCES
---------
-0- DOI reference of this software release
https://zenodo.org/badge/latestdoi/23060862


-1- If the Finite Difference code has helped you in your research please refer to the following paper in your publications:


Finite-difference modeling experiments for seismic interferometry
Jan Thorbecke and Deyan Draganov
2011, Geophysics, Vol. 76, no. 6 (November-December); p H1--H18, doi: 10.1190/GEO2010-0039.1
Download: https://janth.home.xs4all.nl/Publications/Articles/ThorbeckeDraganov2012.pdf



-2- If the Marchenko code has helped you in your research please refer to this paper in your publications:


Implementation of the Marchenko method
Jan Thorbecke, Evert Slob, Joeri Brackenhoff, Joost van der Neut, and Kees Wapenaar
2017, Geophysics, Vol. 82, no. 6 (November-December); p. WB29--WB45, doi: 10.1190/GEO2017-0108.1
Download: https://janth.home.xs4all.nl/Publications/Articles/ThorbeckeGPY2017.pdf


-3- If you used the code to construct homogeneous Green's functions, please refer to this paper in your related publications:


Virtual acoustics in inhomogeneous media with single-sided access:
 Wapenaar, K., Brackenhoff, J., Thorbecke, J., van der Neut, J., Slob, E., and Verschuur, E.,
2018, Scientific Reports, Vol. 8, 2497.
Download: https://www.keeswapenaar.nl/4_Journals/Nature/SR_18.pdf


-4- When you are using the marchenko_primaries algorithm developed by Lele Zhang please refer to the following papers


Free-surface and internal multiple elimination in one step without adaptive subtraction
Lele Zhang and Evert Slob
2019, Geophysics, Vol. 84, no. 1 (January-February); p. A7-A11, doi: 10.1190/GEO2018-0548.1
Download: https://www.keeswapenaar.nl/BeyondInterferometry/geo_19h.pdf


and


Implementation of the Marchenko Multiple Elimination algorithm,
Jan Thorbecke, Lele Zhang, Kees Wapenaar, Evert Slob,
2021, Geophysics, Vol. 86, no. 2 (March-April); p. 1-15, doi: 10.1190/GEO2020-0196.1
Download: https://www.keeswapenaar.nl/4_Journals/Geophysics/geo_21b.pdf


-5- If you use the fdacrtmc code of Max Holicki please refer to the following paper:


Acoustic directional snapshot wavefield decomposition
Holicki, M., Drijkoningen, G., and Wapenaar, K., 2019, Acoustic directional snapshot wavefield decomposition:
Geophysical Prospecting, Vol. 67, 32-51
Download: https://www.keeswapenaar.nl/4_Journals/Geophys.Prosp/GP_19a.pdf


-6- If you use the vmar code of Johno van IJsseldijk please refer to the following paper:


Extracting small time-lapse traveltime changes in a reservoir using primaries and internal multiples after Marchenko-based target zone isolation:
Van IJsseldijk, J., van der Neut, J., Thorbecke, J., and Wapenaar, K., 2023,
Geophysics, Vol. 88 (2), R135-R143.
Download: https://keeswapenaar.nl/4_Journals/Geophysics/geo_23a.pdf


-7- A reference to the extrapolation and migration programs is


Design of one-way wavefield extrapolation operators, using smooth functions in WLSQ optimisation.
Jan Thorbecke, Kees Wapenaar, Gerd Swinnen,
2004, Geophysics Vol 69, p. 1037-1045
Download: https://www.keeswapenaar.nl/4_Journals/Geophysics/geo_04.pdf


-8- The Marchenko plane-wave algorithm is explained in

Design, implementation and application of the Marchenko Plane-wave algorithm.
Jan Thorbecke, Mohammed Almobarak, Johno van IJsseldijk, Joeri Brackenhoff, Giovanni Meles, Kees Wapenaar,
2024 Computers & Geosciences 187 (2024) 105577
Download: https://www.keeswapenaar.nl/4_Journals/Comp.Geosc/CS_24.pdf

The instructions in marchenko/demo/planewave explain how to reproduce the results in this paper.



UPDATES AND LATEST VERSION
--------------------------
The latest version of the source code and manual can be found at:


git clone https://gitlab.com/geophysicsdelft/OpenSource.git


This repository has a mirror on:


https://github.com/JanThorbecke/OpenSource.git



