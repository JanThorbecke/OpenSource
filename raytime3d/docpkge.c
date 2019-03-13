/*                 
 This file is property of the Colorado School of Mines.
 
 Copyright (C) 2007, Colorado School of Mines,
 All rights reserved.
 
 
 Redistribution and use in source and binary forms, with or 
 without modification, are permitted provided that the following 
 conditions are met:
 
 *  Redistributions of source code must retain the above copyright 
 notice, this list of conditions and the following disclaimer.
 *  Redistributions in binary form must reproduce the above 
 copyright notice, this list of conditions and the following 
 disclaimer in the documentation and/or other materials provided 
 with the distribution.
 *  Neither the name of the Colorado School of Mines nor the names of
 its contributors may be used to endorse or promote products 
 derived from this software without specific prior written permission.
 
 Warranty Disclaimer:
 THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
 COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
 STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 POSSIBILITY OF SUCH DAMAGE.
 
 
 Export Restriction Disclaimer:
 We believe that CWP/SU: Seismic Un*x is a low technology product that does
 not appear on the Department of Commerce CCL list of restricted exports.
 Accordingly, we believe that our product meets the qualifications of
 an ECCN (export control classification number) of EAR99 and we believe
 it fits the qualifications of NRR (no restrictions required), and
 is thus not subject to export restrictions of any variety.
 
 Approved Reference Format:
 In publications, please refer to SU as per the following example:
 Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
 Release No. __: an open source software  package for seismic 
 research and processing, 
 Center for Wave Phenomena, Colorado School of Mines.
 
 Articles about SU in peer-reviewed journals:
 Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
 Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
 Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
 Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
 
 Acknowledgements:
 SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
 School of Mines, partially based on Stanford Exploration Project (SEP) 
 software.
 */

/*********************** self documentation **********************/
/*************************************************************************** 
DOCPKGE - Function to implement the CWP self-documentation facility

requestdoc	give selfdoc on user request (i.e. when name of main is typed)
pagedoc		print self documentation string

**************************************************************************** 
Function Prototypes:
void requestdoc(flag);
void pagedoc();

**************************************************************************** 
requestoc:
Input:
flag		integer specifying i.o. cases

pagedoc():
Returns:	the self-documentation, an array of strings

**************************************************************************** 
Notes:
requestdoc:
In the usual case, stdin is used to pass in data.  However,
some programs (eg. synthetic data generators) don't use stdin
to pass in data and some programs require two or more arguments
besides the command itself (eg. sudiff) and don't use stdin.
In this last case, we give selfdoc whenever too few arguments
are given, since these usages violate the usual SU syntax.
In all cases, selfdoc can be requested by giving only the
program name.

The flag argument distinguishes these cases:
            flag = 0; fully defaulted, no stdin
            flag = 1; usual case
            flag = n > 1; no stdin and n extra args required

pagedoc:
Intended to be called by requesdoc(), but conceivably could be
used directly as in:
      if (xargc != 3) selfdoc();

Based on earlier versions by:
SEP: Einar Kjartansson, Stew Levin CWP: Jack Cohen, Shuki Ronen
HRC: Lyle

**************************************************************************** 
Author: Jack K. Cohen, Center for Wave Phenomena
****************************************************************************/
/**************** end self doc ********************************/

#include "par.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif

 
/*  definitions of global variables */
int xargc; char **xargv;


void requestdoc(int flag)
/*************************************************************************** 
print selfdocumentation as directed by the user-specified flag
**************************************************************************** 
Notes:
In the usual case, stdin is used to pass in data.  However,
some programs (eg. synthetic data generators) don't use stdin
to pass in data and some programs require two or more arguments
besides the command itself (eg. sudiff) and don't use stdin.
In this last case, we give selfdoc whenever too few arguments
are given, since these usages violate the usual SU syntax.
In all cases, selfdoc can be requested by giving only the
program name.

The flag argument distinguishes these cases:
            flag = 0; fully defaulted, no stdin
            flag = 1; usual case
            flag = n > 1; no stdin and n extra args required

pagedoc:
Intended to be called by pagedoc(), but conceivably could be
used directly as in:
      if (xargc != 3) selfdoc();

**************************************************************************** 
Authors: Jack Cohen, Center for Wave Phenomena, 1993, based on on earlier
versions by:
SEP: Einar Kjartansson, Stew Levin CWP: Jack Cohen, Shuki Ronen
HRC: Lyle
****************************************************************************/
{
        switch(flag) {
        case 1:
                if (xargc == 1 && isatty(STDIN)) pagedoc();
        break;
        case 0:
                if (xargc == 1 && isatty(STDIN) && isatty(STDOUT)) pagedoc();
        break;
        default:
                if (xargc <= flag) pagedoc();
        break;
        }
        return;
}


void pagedoc(void)
{
        extern char *sdoc[];
	char **p = sdoc;
        FILE *fp;

        fflush(stdout);
        fp = popen("more -22 1>&2", "w");
	while(*p) (void)fprintf(fp, "%s\n", *p++);
        pclose(fp);

        exit(EXIT_FAILURE);
}
/*----------------------End of Package--------------------------------*/
