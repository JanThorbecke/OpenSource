
ACKNOWLEDGEMENT
---------
This work received funding from the European Research Council (grant 742703) and the  NWO Domain Applied and Engineering Sciences (grant 13939).

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

-2- If the Machenko code has helped you in your research please refer to this paper in your publications:

Implementation of the Marchenko method
Jan Thorbecke, Evert Slob, Joeri Brackenhoff, Joost van der Neut, and Kees Wapenaar
2017, Geophysics, Vol. 82, no. 6 (November-December); p. WB29--WB45, doi: 10.1190/GEO2017-0108.1
Download: https://janth.home.xs4all.nl/Publications/Articles/ThorbeckeGPY2017.pdf

-3- If you used the code to construct homogenoeus Green's functions, please refer to this paper in your related publications:

Virtual acoustics in inhomogeneous media with single-sided access:
 Wapenaar, K., Brackenhoff, J., Thorbecke, J., van der Neut, J., Slob, E., and Verschuur, E., 
2018, Scientific Reports, Vol. 8, 2497. 
Download: http://homepage.tudelft.nl/t4n4v/4_Journals/Nature/SR_18.pdf

-4- When you are using the marchenko_primaries algorithm developed by Lele Zhang please refer to the following papers

Free-surface and internal multiple elimination in one step without adaptive subtraction
Lele Zhang and Evert Slob
2019, Geophysics, Vol. 84, no. 1 (January-February); p. A7-A11, doi: 10.1190/GEO2018-0548.1
Download: http://homepage.tudelft.nl/t4n4v/BeyondInterferometry/geo_19h.pdf

and

Implementation of the Marchenko Multiple Elimination algorithm,
Jan Thorbecke, Lele Zhang, Kees Wapenaar, Evert Slob,
2021, Geophysics, Vol. 86, no. 2 (March-April); p. 1-15, doi: 10.1190/GEO2020-0196.1

-5- If you use the fdacrtmc code of Max Holicki please refer to the following paper:

Acoustic directional snapshot wavefield decomposition
Holicki, M., Drijkoningen, G., and Wapenaar, K., 2019, Acoustic directional snapshot wavefield decomposition: Geophysical
Prospecting, Vol. 67, 32-51
Download: http://homepage.tudelft.nl/t4n4v/4_Journals/Geophys.Prosp/GP_19a.pdf


INSTALLATION
-------------
See the seperate INSTALL file.


REPRODUCING
----------
Almost all Figures in the papers mentioned above can be reproduced with the sofwtare in this repository. Please see the file REPRODUCE for further instructions



Finite Difference Modeling: FDELMODC
------------------------------------
If the compilation has finished without errors and produced an executable called bin/fdelmodc you can run one of the demo programs by running

> ./fdelmodc_plane.scr

in the directory fdelmodc/demo/ 

The demo directory contains many scripts which demonstrate the different possibilities of the modeling program. 

To reproduce the Figures shown in the GEOPHYICS manuscript "Finite-difference modeling experiments for seismic interferometry" the scripts in FiguresPaper directory can be used. Please read the README in the FiguresPaper directory for more instructions and guidelines. 

An extensive manual of fdelmodc can be found in doc/fdelmodcManual.pdf

Marchenko method : MARCHENKO
----------------------------
If the compilation has finished without errors and produced an executable called bin/marchenko you can run one of the demo programs by running a set of scripts that are explained in a README in one of the directories marchenko/demo/oneD or demo/twoD

To reproduce the Figures shown in the GEOPHYICS paper "Implementation of the Marchenko method" the scripts in marchenko/demo/oneD directory can be used. The README in this directory gives more instructions and guidelines. 

To reproduce the Figures shown in the Scientific Reports paper "Virtual acoustics in inhomogeneous media with single-sided access" the scripts in marchenko/demo/ScientificReports directory can be used. The README in this directory gives more instructions and guidelines. 

To reproduce the Figures shown in the GEOPHYICS paper "Implementation of the Marchenko Multiple Elimination algorithm" the scripts in marchenko/demo/mme directory can be used. The README_PRIMARIES in this directory gives more instructions and guidelines. 
A bried manual about the MME program 'marchencko_primaries' can be found in doc/MMEmanual.pdf


MDD
---
The MDD kernels depend on BLAS and LAPACK calls. Free downloads of these libraries can be found on 

https://www.netlib.org/blas/index.html
https://www.netlib.org/lapack/index.html


MKL libraries
-------------
If you are running on x64_86 processors you can download (for free) Intel's highly optimised MKL package:

https://software.intel.com/en-us/mkl/choose-download
 
These Libraries include highly optimised libraries of BLAS, LAPACK, and FFT(W). 
Usually the MKL libraries are installed in $MKLROOT. If that variable is not set, you can try to find the correct path by searching for one of the libraries:

find /opt/intel -name libmkl_gf_lp64.so

and adjust MKLROOT in Make_include accordingly. 
You can also completely disable the use of MKL by commenting out the MKL parts in Make_include. 

In case MKL is installed on your system there is no need to install the netlib packages mentioned for MDD. 


SEISMIC UNIX
-----------
If you want to use the .su files with SU from CWP:
git clone https://github.com/JohnWStockwellJr/SeisUnix
git clone https://gitlab.tudelft.nl/Geophysics/OpenSource.git

==> Please make sure that SU is compiled without XDR (in $CWPROOT/Makefile.config make sure that XDRFLAG is NOT enabled). The SU output files of fdelmodc are all based on local IEEE data.
To exclude the XDRFLAG in SU you have to use the following line in $CWPROOT/src/Makefile.config around line 35:

XDRFLAG =

so the XDRFLAG is empty. If you have made that change you have to remake SU with the commands:
make remake
make xtremake


ZFP
---
The fdacrtmc and the 3D Marchenko code makes use of ZFP compression to store snaphots in CPU memory or do effient IO. This package is included in this repository for your convenience. The latest package and detailed explanation can be found on:

https://github.com/LLNL/zfp

and written by Peter Lindstrom

FDACRTMC
--------
fdacrtmc uses FFTW and the wisdom computations are stored on disk for re-usage.  This directory is defined in fdacrtmc.h

#ifndef WISDOMDIR
#define WISDOMDIR "/tmp/fftw/"
#endif

For the code to run properly the directory /tmp/fftw/ must exist. It is also possible to compile the fdacrtmc code that defines
WISDOMDIR. In fdacrtmc/Makefile add:

CFLAGS  += -DWISDOMDIR="/directory/that/exists"

or you can change the name of WISDOMDIR in fdacrtmc.h

The Finite Different based RTM code (FDACRTMC) will not be installed if MKL is not available. This code depends on FFTW and we do want to
include the full FFTW source code package to this. Our aim is to keep the compilation of the code in the GitHub repository
as simple as possible with the smallest number of external libraries.

UPDATES AND LATEST VERSION
--------------------------
The latest version of the source code and manual can be found at:

git clone https://github.com/JanThorbecke/OpenSource.git

git clone git://github.com/JanThorbecke/OpenSource.git

The code is used by many different people and if there is a request for a new option in the code, then we will try to implement, test and make it available. 

