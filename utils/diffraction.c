#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  insert diffractor in the model used, in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void makesquare(int ix, int iz, int diffrwidth, float **med, float value)
{
    int j,k;
    for (j=-diffrwidth/2; j<diffrwidth/2; j++){
        for (k=-diffrwidth/2; k<diffrwidth/2; k++){
            med[ix+j][iz+k] = value;
        }
    }
    return;
}

void makediamond(int ix, int iz, int diffrwidth, float **med, float value)
{
    int j, k, id;
    
    id=diffrwidth/2;
    j=-id; k=0;
    med[ix+j][iz+k] = value;
    j=+id; k=0;
    med[ix+j][iz+k] = value;
    for (j=-id+1; j<=0; j++){
        for (k=-(j+id); k<=j+id; k++){
            med[ix+j][iz+k] = value;
        }
    }
    for (j=1; j<=id-1; j++){
        for (k=j-id; k<=-(j-id); k++){
            med[ix+j][iz+k] = value;
        }
    }

    return;
}
    
void makecircle(int ix, int iz, int diffrwidth, float **med, float value)
{
    int j,k, id, ir;
    
    id=diffrwidth/2;
    for (j=-id; j<=id; j++){
        for (k=-id; k<=id; k++){
            ir = NINT(sqrt(j*j+k*k));
            if (ir <= id) med[ix+j][iz+k] = value;
        }
    }
    return;
}

    
void diffraction(float *x, float *z, int nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, int *zp, int nx, int diffrwidth, int type)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		interface[i] = 0.0;
		zp[i] = 0;
	}

	if (gridcs == NULL && gridro == NULL) {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
            switch( type )
            {
                case 1:
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    break;
                case 2:
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    break;
                default:
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    break;
            }
		}
	}
	else if (gridcs == NULL) {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
            switch( type )
            {
                case 1:
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    break;
                case 2:
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    break;
                default :
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    break;
            }
		}
	}
	else {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
            switch( type )
            {
                case 1:
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    makediamond(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcs, cs[2][i]);
                    break;
                case 2:
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    makecircle(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcs, cs[2][i]);
                    break;
                default :
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcp, cp[2][i]);
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridro, ro[2][i]);
                    makesquare(NINT(x[i]/dx), NINT(z[i]/dz), diffrwidth,  gridcs, cs[2][i]);
                    break;
            }
		}
	}

	return;
}
