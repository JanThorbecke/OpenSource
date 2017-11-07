//
//  JespersRayTracer.c
//  
//
//  Written by Jesper Spetzler
//
//  changed to C by Jan Thorbecke on 21/09/2017.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "raytime.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

static float H, L, W;

typedef struct _icoord { /* 3D coordinate integer */
    int z;
    int x;
    int y;
} icoord;

typedef struct _fcoord { /* 3D coordinate float */
    float z;
    float x;
    float y;
} fcoord;

int getnRay(icoord size, fcoord s, fcoord r, float dx, int nRayStep);
int traceTwoPoint(fcoord s, fcoord r, int nRay, fcoord *rayReference3D);
float takeOffAngle(fcoord s, fcoord r);
float referenceSlowness(float *slowness, icoord size, int nRay, fcoord r, fcoord s);
int xPointIndex(const float _x, int nx, float L);
int zPointIndex(const float _z, int nz, float H);
int yPointIndex(const float _y, int ny, float W);
fcoord getSlownessGradient(const float _x, const float _z, float *slowness, icoord size);
float qMulGradU1(const float _x, const float _z, const float _angle, float *slowness, icoord size);
float greenTwoP(const float _so, const float _slow, const float _sL, int nRay, fcoord s, fcoord r, float *slowness, icoord size);
float qatso(const float _so, const float _angle, int nRay, fcoord s, fcoord r, fcoord *rayReference3D, float *slowness, icoord size);
float slownessA(float *slowness, icoord size, float _x, float _y, float _z);
float getdT2(const float _x, const float _z, const float so, const float _angle, const float _ds, int nRay, fcoord s, fcoord r, fcoord *rayReference3D, float *slowness, icoord size);
float greenIntP(const float _so, const float _s, const float _sL, float *slowness, icoord size, int nRay, fcoord r, fcoord s);
float secondDerivativeU1(float *slowness, icoord size, const float _x, const float _z, const float _angle, fcoord s, fcoord r);
int calculatePerturbedRay(fcoord *rayPerturbed3D, fcoord s, fcoord r, int nRay, fcoord *rayReference3D, float *slowness, icoord size);
float angle2qx(const float _angle);
float angle2qz(const float _angle);
float ModelInterpolation_slowness2D(float *slowness, icoord size, const float _x, const float _z);
float ModelInterpolation_slowness3D(float *slowness, icoord size, const float _x, const float _z, const float _y);
void applyMovingAverageFilter(float *slowness, icoord size, int window, int dim, float *averageModel);



#define lGradient 1
#define EPSMIN 0.1
#define minValueGradient 1e-10
#define PI 3.1514926535
#define minValueSecondDerivativeU1 1e-6
#define DPHI_ANGLE 1.0   // 0.5

int getWaveParameter(float *slowness, icoord size, float dgrid, fcoord s, fcoord r, rayPar ray, fcoord *T, float *Jr)
{
	static int first=1;
	float *smooth;
    float T0, T1, T2;
    float uo, u1, lengthRefRay;
    float x, y, z;
    float dx, dy, dz, dl, so, ds;
    float angle;
    float dQdPhi, J, greentmp;
    int nRayTmp, error, i;
    fcoord *rayReference3D;
    
    T0 = T1 = T2 = 0;
    J = 1;
    error = 0;
    
    nRayTmp = getnRay(size, s, r, dgrid, ray.nray);
    
    //fprintf(stderr,"Calling getnRay gives nRayTmp=%d nRayStep=%d\n", nRayTmp, nRayStep);

    rayReference3D = (fcoord *)calloc(nRayTmp,sizeof(fcoord));
    traceTwoPoint(s, r, nRayTmp, rayReference3D);
    
    dx = rayReference3D[nRayTmp-1].x - rayReference3D[0].x;
    dy = rayReference3D[nRayTmp-1].y - rayReference3D[0].y;
    dz = rayReference3D[nRayTmp-1].z - rayReference3D[0].z;
    lengthRefRay = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    
    angle = takeOffAngle(s, r);
    
    if ((lengthRefRay <= 0) || (nRayTmp <= 1))
        return(-1);
    
    uo = referenceSlowness(slowness, size, nRayTmp, r, s);
    
    T0 = lengthRefRay*uo;
    ds = lengthRefRay/(nRayTmp-1);
    J = lengthRefRay;
    dQdPhi = 0;

    for (i = 0; i < nRayTmp-1; i++)
    {
        x = 0.5*(rayReference3D[i+1].x + rayReference3D[i].x);
        y = 0.5*(rayReference3D[i+1].y + rayReference3D[i].y);
        z = 0.5*(rayReference3D[i+1].z + rayReference3D[i].z);
        
        u1 = slownessA(slowness, size, x, z, y) -  uo;
        
        dx = rayReference3D[i+1].x - rayReference3D[i].x;
        dy = rayReference3D[i+1].y - rayReference3D[i].y;
        dz = rayReference3D[i+1].z - rayReference3D[i].z;
        
        dl = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        
        T1 += dl*u1;
        
        so = i*ds;
        
        if (ray.useT2 != 0)
            T2 += getdT2(x, z, so, angle, ds, nRayTmp, s, r, rayReference3D, slowness, size);

        if (ray.geomspread != 0) {
            if (so <= 0) {
                dQdPhi = 0;
            }
            else {
                greentmp = greenIntP(lengthRefRay, so, lengthRefRay, slowness, size, nRayTmp, r, s);
                dQdPhi += greentmp*secondDerivativeU1(slowness, size, x, z, angle, r, s)*ds/so;
            }
        }
    }

    if (ray.useT2)
        T2 *= 0.5;
    
    T->x = T0;
    T->y = T1;
    T->z = T2;
    
    // The geometrical spreading factor
    
    if (ray.geomspread)
    {
        J += dQdPhi;
        
        if (J == 0)
            J = 1;
        
        if (J < 0)
        {
            error = -1; //snegativeGeometricalSpreadingFactor;
            J = fabs(J);
        }
    }
    
    if (size.y == 1) {
        J = sqrt(J);
    }
    
    *Jr = J;
    free(rayReference3D);
    
    return(error);
}

int getnRay(icoord size, fcoord s, fcoord r, float dx, int nRayStep)
{
    int dn, nRayTmp;
    float dl, dr;
    
    H = (size.z-1)*dx;
    L = (size.x-1)*dx;
    W = (size.y-1)*dx;
    
    if (size.y == 1) { // 2D model
        dn = (size.x + size.z)/2;
        dl = sqrt(pow(L, 2) + pow(H, 2))/dn;
        dr = sqrt(pow(r.x-s.x, 2) + pow(r.z-s.z, 2));
    }
    else { // 3D model
        dn = (size.x + size.z + size.y)/3;
        dl = sqrt(pow(L, 2) + pow(H, 2) + pow(W, 2))/dn;
        dr = sqrt(pow(r.x-s.x, 2) + pow(r.z-s.z, 2) + pow(r.y-s.y, 2));
        
    }
    nRayTmp = MIN(300,dr*nRayStep/dl);
    //fprintf(stderr,"getnRay: gives nRayTmp=%d dr=%f dl=%f\n", nRayTmp, dr, dl);

    if (nRayTmp <= nRayStep)
        nRayTmp = nRayStep;
    
    return nRayTmp;

}

int traceTwoPoint(fcoord s, fcoord r, int nRay, fcoord *rayReference3D)
{
    float x, y, z;
    int i;
    
    for (i = 0; i < nRay; i++)
    {
        x = s.x + (r.x - s.x)*i/(nRay-1);
        y = s.y + (r.y - s.y)*i/(nRay-1);
        z = s.z + (r.z - s.z)*i/(nRay-1);
        rayReference3D[i].z=z;
        rayReference3D[i].x=x;
        rayReference3D[i].y=y;
    }
    
    return 0;
}


int calculatePerturbedRay(fcoord *rayPerturbed3D, fcoord s, fcoord r, int nRay, fcoord *rayReference3D, float *slowness, icoord size)
{
    float si, sl, deltaS, gso, angle, qx, qz;
    int i;
    
    sl = sqrt(pow((r.x-s.x), 2) + pow((r.y-s.y), 2) + pow((r.z-s.z), 2));
    
    if ((sl <= 0) || (nRay <= 1))
        return 0;
    
    deltaS = sl/(nRay-1);
    angle = takeOffAngle(s, r);
    
    qx = angle2qx(angle);
    qz = angle2qz(angle);
    
    for (i = 0; i < nRay; i++)
    {
        si = i*deltaS;
        
        gso = qatso(si, angle, nRay, s, r, rayReference3D, slowness, size);
        
        rayPerturbed3D[i].x = rayReference3D[i].x + qx*gso;
        rayPerturbed3D[i].z = rayReference3D[i].z + qz*gso;
        rayPerturbed3D[i].y = rayReference3D[i].y;

    }
    
    return 0;
}

float takeOffAngle(fcoord s, fcoord r)
{
    float angle = 0;
    
    if ((s.x == r.x) && (s.z == r.z))
        angle = PI/2;
    else if ((s.x <= r.x) && (s.z < r.z))
        angle = atan(fabs(r.x-s.x)/fabs(r.z-s.z));
    else if ((s.x < r.x) && (s.z >= r.z))
        angle = PI/2 + atan(fabs(r.z-s.z)/fabs(r.x-s.x));
    else if ((s.x >= r.x) && (s.z > r.z))
        angle = PI + atan(fabs(r.x-s.x)/fabs(r.z-s.z));
    else if ((s.x > r.x) && (s.z <= r.z))
        angle = 3*PI/2 + atan(fabs(r.z-s.z)/fabs(r.x-s.x));
    
    return (angle);
}

float angle2qx(const float _angle)
{
    float qx = 0;
    
    if ((_angle >= 0) && (_angle < PI/2))
        qx = -cos(_angle);
    else if ((_angle >= PI/2) && (_angle < PI))
        qx = sin(_angle - PI/2);
    else if ((_angle >= PI) && (_angle < 3*PI/2))
        qx = cos(_angle - PI);
    else if ((_angle >= 3*PI/2) && (_angle <= 2*PI))
        qx = -sin(_angle - 3*PI/2);
    
    return (qx);
}

float angle2qz(const float _angle)
{
    float qz = 0;
    
    if ((_angle >= 0) && (_angle < PI/2))
        qz = sin(_angle);
    else if ((_angle >= PI/2) && (_angle < PI))
        qz = cos(_angle - PI/2);
    else if ((_angle >= PI) && (_angle < 3*PI/2))
        qz = -sin(_angle - PI);
    else if ((_angle >= 3*PI/2) && (_angle <= 2*PI))
        qz = -cos(_angle - 3*PI/2);
    
    return (qz);
}

// Sofar used in 2D only

float qatso(const float _so, const float _angle, int nRay, fcoord s, fcoord r, fcoord *rayReference3D, float *slowness, icoord size)
{
    float slow, sl, deltaS, x, z;
    float qatsol;
    int i;
    
    sl = sqrt(pow((r.x-s.x),2) + pow((r.z-s.z),2) + pow((r.y-s.y),2));
    
    if ((sl <= 0) || (nRay <= 1))
    {
        return 0;
    }
    
    deltaS = sl/(nRay-1);
    
    qatsol = 0;
    for (i = 0; i < nRay; i++)
    {
        slow = i*deltaS;
        x = rayReference3D[i].x;
        z = rayReference3D[i].z;
//        fprintf(stderr,"qatso: calling greenTwoP for iray %d (/%d)\n",i,nRay);

        qatsol += greenTwoP(_so, slow, sl, nRay, s, r, slowness, size)*qMulGradU1(x, z, _angle, slowness, size)*deltaS;
    }

    return(qatsol);
}

float getdT2(const float _x, const float _z, const float _so, const float _angle, const float _ds, int nRay, fcoord s, fcoord r, fcoord *rayReference3D, float *slowness, icoord size)
{
    float T2 = 0;
    float qatsol;
    float qMulGradU1l;

 //   fprintf(stderr,"getdT2: calling qatso nRay=%d\n",nRay);

    qatsol = qatso(_so, _angle, nRay, s, r, rayReference3D, slowness, size);
    
//    fprintf(stderr,"getdT2: calling qMulGradU1\n");

    qMulGradU1l = qMulGradU1(_x, _z, _angle, slowness, size);

    T2 = qatsol*qMulGradU1l*_ds;
    
    return(T2);
}

float greenTwoP(const float _so, const float _slow, const float _sL, int nRay, fcoord s, fcoord r, float *slowness, icoord size)
{
    float greenTwoP = 0;
    float uo = referenceSlowness(slowness, size, nRay, r, s);
    
//    fprintf(stderr,"greenTwoP: slowness = %f nRay=%d\n",uo,nRay);

    if (_sL <= 0)
    {
        return(0);
    }
    
    if (_slow <= _so)
        greenTwoP = -(1 - _so/_sL)*_slow/uo;
    else
        greenTwoP = -_so*(1-_slow/_sL)/uo;
    
    return(greenTwoP);
}

float qMulGradU1(const float _x, const float _z, const float _angle, float *slowness, icoord size)
{
    float qMulGradU1;
    float gradu1x, gradu1z;
    float qx, qz;
    fcoord slownessGradient;
    
    slownessGradient = getSlownessGradient(_x, _z, slowness, size);
    gradu1x = slownessGradient.x;
    gradu1z = slownessGradient.z;
    
    qx = angle2qx(_angle);
    qz = angle2qz(_angle);
    
    qMulGradU1 = qx*gradu1x + qz*gradu1z;
    
    return(qMulGradU1);
}

float referenceSlowness(float *slowness, icoord size, int nRay, fcoord r, fcoord s)
{
    float x, y, z;
    float uo = 0;
    int i;
    
    for (i = 0; i < nRay; i++)
    {
        x = s.x + (r.x - s.x)*i/(nRay-1);
        z = s.z + (r.z - s.z)*i/(nRay-1);

        if (size.y == 1) // 2D
            uo += ModelInterpolation_slowness2D(slowness, size, x, z);
        else
        {
            y = s.y + (r.y - s.y)*i/(nRay-1);
            uo += ModelInterpolation_slowness3D(slowness, size, x, z, y);
        }
    }
    
    uo /= nRay;
    
    return(uo);
}

fcoord getSlownessGradient(const float _x, const float _z, float *slowness, icoord size)
{
    float dx, dz, x1, x2, z1, z2;
    float slow2, slow1;
    float gradu1x, gradu1z;
    fcoord slownessGradient;
    
    dx = lGradient*L/(size.x-1);
    dz = lGradient*H/(size.z-1);
    
    x1 = _x-dx;
    x2 = _x+dx;
    
    if (x1 <= 0)
        x1 = EPSMIN;
    
    if (x2 >= L)
        x2 = L - EPSMIN;
    
    if (size.y == 1)
    {
        slow1 = ModelInterpolation_slowness2D(slowness, size, x1, _z);
        slow2 = ModelInterpolation_slowness2D(slowness, size, x2, _z);
    }
    else
    {
        slow1 = ModelInterpolation_slowness3D(slowness, size, x1, _z, 0);
        slow2 = ModelInterpolation_slowness3D(slowness, size, x2, _z, 0);
    }
    
    if (fabs(slow2-slow1) < minValueGradient)
        gradu1x = 0;
    else
        gradu1x = (slow2 - slow1)/(x2-x1);
    
    z1 = _z-dz;
    z2 = _z+dz;
    
    if (z1 <= 0)
        z1 = EPSMIN;
    
    if (z2 >= H)
        z2 = H - EPSMIN;
    
    if (size.y == 1)
    {
        slow1 = ModelInterpolation_slowness2D(slowness, size, _x, z1);
        slow2 = ModelInterpolation_slowness2D(slowness, size, _x, z2);
    }
    else
    {
        slow1 = ModelInterpolation_slowness3D(slowness, size, _x, z1, 0);
        slow2 = ModelInterpolation_slowness3D(slowness, size, _x, z2, 0);
    }
    
    if (fabs(slow2-slow1) < minValueGradient)
        gradu1z = 0;
    else
        gradu1z = (slow2 - slow1)/(z2-z1);
    
    slownessGradient.x=gradu1x;
    slownessGradient.z=gradu1z;
    slownessGradient.y=0;

    return(slownessGradient);
}

int xPointIndex(const float _x, int nx, float L)
{
    int i;
    
    if (_x <= 0)
        return(0);
    
    if (_x >= L)
        i = nx - 1;
    else
    {
        if (0 < L)
            i = _x*nx/L;
        else
            i = 0;
    }
    
    return(i);
}

int zPointIndex(const float _z, int nz, float H)
{
    int i;
    
    if (_z <= 0) return(0);
    
    if (_z >= H)
        i = nz - 1;
    else
    {
        if (0 < H)
            i = _z*nz/H;
        else
            i = 0;
    }
    
    return(i);
}

int yPointIndex(const float _y, int ny, float W)
{
    int i;
    
    if (_y <= -0.5*W)
        return(0);
    
    if (_y >= 0.5*W)
        i = ny - 1;
    else
    {
        if (0 < W)
            i = ny*(_y/W + 0.5);
        else
            i = 0;
    }
    
    return(i);
}

float ModelInterpolation_slowness2D(float *slowness, icoord size, const float _x, const float _z)
{
    float slow;
    float f11, f12, f21, f22;
    float t, j;
    float x1, x2;
    float z1, z2;
    int nx, nz, ix, iz, ixMin, ixMax, izMin, izMax;
    int ixCoordinate, izCoordinate;
    
    slow = f11 = f12 = f21 = f22 = 0;
    nx = size.x;
    nz = size.z;
    
    ixCoordinate = (int) _x*nx/L;

    if (ixCoordinate >= nx)
        ixCoordinate = nx;
    
    if (ixCoordinate == nx)
    {
        x1 = (float) L*(ixCoordinate-1)/nx;
        x2 = (float) L;
    }
    else if (ixCoordinate <= 0)
    {
        x1 = 0;
        x2 = (float) L/nx;
    }
    else
    {
        x1 = (float) L*ixCoordinate/nx;
        x2 = (float) L*(ixCoordinate+1)/nx;
    }
    
    if (x1 < 0)
        x1 = 0;
    
    if (x1 > L)
        x1 = L;
    
    if (x2 < 0)
        x2 = 0;
    
    if (x2 > L)
        x2 = L;
    
    izCoordinate = (int) _z*nz/H;
    
    if (izCoordinate >= nz)
        izCoordinate = nz;
    
    if (izCoordinate == nz)
    {
        z1 = (float) H*(izCoordinate-1)/nz;
        z2 = (float) H;
    }
    else if (izCoordinate <= 0)
    {
        z1 = 0;
        z2 = (float) H/nz;
    }
    else
    {
        z1 = (float) H*izCoordinate/nz;
        z2 = (float) H*(izCoordinate+1)/nz;
    }
    
    if (z1 < 0)
        z1 = 0;
    
    if (z1 > H)
        z1 = H;
    
    if (z2 < 0)
        z2 = 0;
    
    if (z2 > H)
        z2 = H;
    
    ix = xPointIndex(_x, size.x, L);
    iz = zPointIndex(_z, size.z, H);
    
    if (ix == 0)
    {
        ixMin = 0;
        ixMax = 1;
    }
    else if (ix == nx-1)
    {
        ixMin = nx-2;
        ixMax = nx-1;
    }
    else
    {
        ixMin = ix-1;
        ixMax = ix+1;
    }
    
    if (iz == 0)
    {
        izMin = 0;
        izMax = 1;
    }
    else if (iz == nz-1)
    {
        izMin = nz-2;
        izMax = nz-1;
    }
    else
    {
        izMin = iz-1;
        izMax = iz+1;
    }
    
    f11 = slowness[ixMin*size.z+izMin];
    f21 = slowness[ixMax*size.z+izMin];
    f12 = slowness[ixMin*size.z+izMax];
    f22 = slowness[ixMax*size.z+izMax];
    
    t = (_x-x1)/(x2-x1);
    j = (_z-z1)/(z2-z1);
    
    slow = f11*(1-t)*(1-j) + f21*t*(1-j) + f12*(1-t)*j + f22*t*j;
    
    return (slow);
}

float ModelInterpolation_slowness3D(float *slowness, icoord size, const float _x, const float _z, const float _y)
{
    float slow;
    float f111, f112, f212, f211;
    float f121, f122, f222, f221;
    float t, j, r;
    float x1, x2;
    float y1, y2;
    float z1, z2;
    int ix, iy, iz, ixMin, ixMax, iyMin, iyMax, izMin, izMax;
    int nx, nz, ny, nxz;
    int ixCoordinate, iyCoordinate, izCoordinate;
    
    nx = size.x;
    nz = size.z;
    ny = size.y;
    nxz = nx*nz;

    slow = f111 = f112 = f212 = f211 = f121 = f122 = f222 = f221 = 0;
    
    ixCoordinate = _x*nx/L;
    
    if (ixCoordinate >= nx)
        ixCoordinate =  nx;
    
    if (ixCoordinate == nx)
    {
        x1 = (float) L*(ixCoordinate-1)/nx;
        x2 = L;
    }
    else if (ixCoordinate <= 0)
    {
        x1 = 0;
        x2 = (float) L/nx;
    }
    else
    {
        x1 = (float) L*ixCoordinate/nx;
        x2 = (float) L*(ixCoordinate+1)/nx;
    }
    
    if (x1 < 0)
        x1 = 0;
    
    if (x1 > L)
        x1 = L;
    
    if (x2 < 0)
        x2 = 0;
    
    if (x2 > L)
        x2 = L;
    
    izCoordinate = _z*nz/H;
    
    if (izCoordinate >= nz)
        izCoordinate = nz;
    
    if (izCoordinate == nz)
    {
        z1 = H*(izCoordinate-1)/nz;
        z2 = H;
    }
    else if (izCoordinate <= 0)
    {
        z1 = 0;
        z2 = (float) H/nz;
    }
    else
    {
        z1 = (float) H*izCoordinate/nz;
        z2 = (float) H*(izCoordinate+1)/nz;
    }
    
    if (z1 < 0)
        z1 = 0;
    
    if (z1 > H)
        z1 = H;
    
    if (z2 < 0)
        z2 = 0;
    
    if (z2 > H)
        z2 = H;
    
    iyCoordinate = ny*(_y/W + 0.5);
    
    if (iyCoordinate >= ny)
        iyCoordinate = ny;
    
    if (iyCoordinate == ny)
    {
        y1 = (float) W*(iyCoordinate-1-0.5*ny)/ny;
        y2 = 0.5*W;
    }
    else if (iyCoordinate <= 0)
    {
        y1 = -0.5*W;
        y2 = (float) W*(1-0.5*ny)/ny;
    }
    else
    {
        y1 = (float) W*(iyCoordinate-0.5*ny)/ny;
        y2 = (float) W*(iyCoordinate+1-0.5*ny)/ny;
    }
    
    if (y1 < -0.5*W)
        y1 = -0.5*W;
    
    if (y1 > 0.5*W)
        y1 = 0.5*W;
    
    if (y2 < -0.5*W)
        y2 = -0.5*W;
    
    if (y2 > 0.5*W)
        y2 = 0.5*W;
    
    ix = xPointIndex(_x, size.x, L);
    iy = yPointIndex(_y, size.y, W);
    iz = zPointIndex(_z, size.z, H);
    
    if (ix == 0)
    {
        ixMin = 0;
        ixMax = 1;
    }
    else if (ix == nx-1)
    {
        ixMin = nx-2;
        ixMax = nx-1;
    }
    else
    {
        ixMin = ix-1;
        ixMax = ix+1;
    }
    
    if (iz == 0)
    {
        izMin = 0;
        izMax = 1;
    }
    else if (iz == nz-1)
    {
        izMin = nz-2;
        izMax = nz-1;
    }
    else
    {
        izMin = iz-1;
        izMax = iz+1;
    }
    
    if (iy == 0)
    {
        iyMin = 0;
        iyMax = 1;
    }
    else if (iy == ny-1)
    {
        iyMin = ny-2;
        iyMax = ny-1;
    }
    else
    {
        iyMin = iy-1;
        iyMax = iy+1;
    }
    
    nxz = nx*nz;
    f111 = slowness[iyMin*nxz+ixMin*nz+izMin];
    f211 = slowness[iyMax*nxz+ixMin*nz+izMin];
    f121 = slowness[iyMin*nxz+ixMax*nz+izMin];
    f221 = slowness[iyMax*nxz+ixMax*nz+izMin];
    f112 = slowness[iyMin*nxz+ixMin*nz+izMax];
    f212 = slowness[iyMax*nxz+ixMin*nz+izMax];
    f122 = slowness[iyMin*nxz+ixMax*nz+izMax];
    f222 = slowness[iyMax*nxz+ixMax*nz+izMax];
    
    //    cout << "slowness3D 6 "  << endl;
    
    r = (_z-z1)/(z2-z1);
    t = (_x-x1)/(x2-x1);
    j = (_y-y1)/(y2-y1);
    
    slow = f111*(1-t)*(1-j)*(1-r) + f112*(1-t)*(1-j)*r + f211*t*(1-j)*(1-r) + f212*t*(1-j)*r + f121*(1-t)*j*(1-r) + f122*(1-t)*j*r + f222*t*j*r + f221*t*j*(1-r);
    
    slow = f111*(1-r)*(1-t)*(1-j) + f112*(1-r)*(1-t)*j + f211*r*(1-t)*(1-j) + f212*r*(1-t)*j + f121*(1-r)*t*(1-j) + f122*(1-r)*t*j + f222*r*t*j + f221*r*t*(1-j);
    
    
    //    if (slow != slow)
    /*
    if (slow <= 0)
    {
        cout << "                  ModelInterpolation::slowness3D " << 1/slow << "  " << 1/f111 << "  " << 1/f112 << "  " << 1/f211 << "  " << 1/f212 << "  " << 1/f121 << "  " << 1/f122 << "  " << 1/f222 << "  " << 1/f211 << "  " << r << "  " << t << "  " << j << "  " << ixCoordinate << "  "  << x1 << "  " << x2 << "  " << _x << "  " << nx << "  "  << L <<  endl;
        cout << "                        ModelInterpolation::slowness3D, x1, x2 = " << x1 << "  " << x2 << endl;
        cout << "                        ModelInterpolation::slowness3D, y1, y2 = " << y1 << "  " << y2 << endl;
        cout << "                        ModelInterpolation::slowness3D, z1, z2 = " << z1 << "  " << z2 << "  " << _z << endl;
        
        
        exit(EXIT_FAILURE);
    }
    */
    
    return (slow);
}

float slownessA(float *slowness, icoord size, float _x, float _z, float _y)
{
    float slow;
    
    if (size.y == 1)
        slow = ModelInterpolation_slowness2D(slowness, size, _x, _z);
    else
        slow = ModelInterpolation_slowness3D(slowness, size, _x, _z, _y);
    
    return(slow);
}

float greenIntP(const float _so, const float _s, const float _sL, float *slowness, icoord size, int nRay, fcoord r, fcoord s)
{
    float greenIntP;
    float uo = referenceSlowness(slowness, size, nRay, r, s);
    
    if (_sL <= 0)
    {
        greenIntP = 0;
        return(greenIntP);
    }
    
    if (_s <= _so)
        greenIntP = (_so - _s)/uo;
    else
        greenIntP = 0;
    
    return(greenIntP);
}

float secondDerivativeU1(float *slowness, icoord size, const float _x, const float _z, const float _angle, fcoord r, fcoord s)
{
    float secondDerivativeU1 = 0;
    float dphi, sl;
    float qx, qz;
    float dh, x1, z1, x2, z2;
    
    dphi = DPHI_ANGLE*PI/180.0;
    sl = sqrt(pow((r.x-s.x),2) + pow((r.z-s.z),2) + pow((r.y-s.y),2));
    
    // Here qx and qz are perpendicular to the raz direction
    
    qx = angle2qx(_angle);
    qz = angle2qz(_angle);
    
    dh = sl*tan(2*dphi);
    x2 = _x + dh*qx;
    z2 = _z + dh*qz;
    
    x1 = _x - dh*qx;
    z1 = _z - dh*qz;
    
    if (x1 <= 0)
        x1 = EPSMIN;
    
    if (x1 >= L)
        x1 = L - EPSMIN;
    
    if (x2 <= 0)
        x2 = EPSMIN;
    
    if (x2 >= L)
        x2 = L - EPSMIN;
    
    if (z1 <= 0)
        z1 = EPSMIN;
    
    if (z1 >= H)
        z1 = H - EPSMIN;
    
    if (z2 <= 0)
        z2 = EPSMIN;
    
    if (z2 >= H)
        z2 = H - EPSMIN;
    
    secondDerivativeU1 = (slownessA(slowness, size, x2, z2, 0) + slownessA(slowness, size, x1, z1, 0) - 2*slownessA(slowness, size, _x, _z, 0))/(4*pow(dphi, 2));
    
    if (fabs(secondDerivativeU1) <= minValueSecondDerivativeU1)
        secondDerivativeU1 = 0;
    
    return(secondDerivativeU1);
}

// Moving average filter
void applyMovingAverageFilter(float *slowness, icoord size, int window, int dim, float *averageModel)
{
    float   averageFilter;
    int     nsamp, iAverageX, iAverageY, iAverageZ;
    int     jWindowX, jWindowZ, jWindowY, ix, iy, iz, nw;

    nw = window;
	if (dim == 2) {
        for (ix = 0; ix < size.x; ix++) {
            for (iz = 0; iz < size.z; iz++) {
                averageFilter =  0;
                nsamp = 0;
                for (jWindowX = -nw; jWindowX <= nw; jWindowX++) {
                    iAverageX = nw + ix + jWindowX;	

//					if (iAverageX < 0) iAverageX = 0;
//                    if (iAverageX > size.x-1) iAverageX = size.x-1;

                    for (jWindowZ = -nw; jWindowZ <= nw; jWindowZ++) {
                        iAverageZ = nw + iz + jWindowZ;

//                        if (iAverageZ < 0) iAverageZ = 0;
//                        if (iAverageZ > size.z-1) iAverageZ = size.z-1;

                        averageFilter += slowness[iAverageX*size.z+iAverageZ];
                        nsamp += 1;
                    }
                }
				if (nsamp > 0) {
                    averageFilter /= nsamp;
				    averageModel[ix*size.z+iz] = averageFilter;
                }
                else 
				    averageModel[ix*size.z+iz] = slowness[(ix+nw)*size.z+iz+nw];
			}
		}
	}
/* 3D ray-tracer not yet implemented 
	else {
	    for (iz = 0; iz < size.z; iz++) {
			for (ix = 0; ix < size.x; ix++) {
                for (iy = 0; iy < size.y; iy++) {
                    averageFilter =  0;
                    nsamp = 0;

					for (jWindowZ = -window; jWindowZ <= window; jWindowZ++) {
                        iAverageZ = iz + jWindowZ;

						if (iAverageZ < 0) iAverageZ = 0;
	                    if (iAverageZ > size.z-1) iAverageZ = size.z-1;

                        for (jWindowX = -window; jWindowX <= window; jWindowX++) {
                            iAverageX = ix + jWindowX;	

							if (iAverageX < 0) iAverageX = 0;
                            if (iAverageX > size.x-1) iAverageX = size.x-1;

                            for (jWindowY = -window; jWindowY <= window; jWindowY++) {
                                iAverageY = iy + jWindowY;

                                if (iAverageY < 0) iAverageY = 0;
                                if (iAverageY > size.y-1) iAverageY = size.y-1;

                                averageFilter += slowness[iAverageZ+iAverageX*size.z+iAverageY*size.z*size.x];
                                nsamp += 1;
                            }
                        }
					}
					
					if (nsamp > 0) {
                        averageFilter /= nsamp;
				        averageModel[iz+ix*size.z+iy*size.z*size.x] = averageFilter;
                    }
                    else {
                       averageModel[iz+ix*size.z+iy*size.z*size.x] = slowness[iz+ix*size.z+iy*size.z*size.x];
                    }
				}
			}
		}
	}	
*/

	return;
}

