#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include "segy.h"

void name_ext(char *filename, char *extension);

char *sdoc[] = {
" readlog - convert ASCII file with depth velocity density columns to SU ",
" ",
" Usage: readlog src_txt=mod272log.txt dz=2 file_base=SA.su ",
" ",
"   src_txt= ................. input ASCII text file with 3 columns.",
"   file_base= ............... base name of the output file(s).",
"   dz= ...................... grid distance in z in meters.",
" ",
" src_txt format is given by depth[m] velocity[m/s] density[10^3 kg/m^3].",
" For example:",
" 0	3076.4	2.01 ",
" 25.9	3076.4	2.01 ",
" 133.8	3076.4	2.01 ",
" 157	3081.5	2.15 ",
" ",
" Basic program not very well tested.",
NULL};

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int main(int argc, char **argv)
{
    char *src_txt, *file_base, filename[256];
    FILE *fp, *fpsu;
    int i, iz, is, n1, nsrc;
    float *dep, *vel, *rho, dz, *trrho, *trvel;
    segy hdr;
    size_t  nwrite;


    initargs(argc,argv);
    requestdoc(0);

    if (!getparstring("src_txt",&src_txt)) src_txt=NULL;
    if(getparstring("file_base",&file_base))

    if (!getparfloat("dz",&dz)) dz=4.0;

    fprintf(stderr,"src_txt=%s dz=%f\n",src_txt, dz);
    if (src_txt!=NULL) {
        /* Open text file */
        nsrc=0;
        fp=fopen(src_txt,"r");
        assert(fp!=NULL);
        /* Get number of lines */
        while (!feof(fp)) if (fgetc(fp)=='\n') nsrc++;
        fseek(fp,-1,SEEK_CUR);
        if (fgetc(fp)!='\n') nsrc++; /* Checks if last line terminated by /n */
        vmess("Number of grid values in src_txt file: %d",nsrc);
        rewind(fp);
    }


    dep = (float *)malloc(nsrc*sizeof(float));
    vel = (float *)malloc(nsrc*sizeof(float));
    rho = (float *)malloc(nsrc*sizeof(float));

    if (src_txt!=NULL) {
        /* Read in source coordinates */
        for (i=0;i<nsrc;i++) {
            if (fscanf(fp,"%e %e %e\n",&dep[i],&vel[i],&rho[i])!=3) vmess("Source Text File: Can not parse coordinates on line %d.",i);
            fprintf(stderr,"Array from src_txt: at %d depth=%f vel=%f rho=%f \n", i, dep[i], vel[i], rho[i]);
        }
        /* Close file */
        fclose(fp);
    }

    n1 = NINT(dep[nsrc-1]/dz)+1;
    vmess("Number of sources in src_txt file: %d with dz=%f Nz=%d",nsrc, dz, n1);

    memset(&hdr,0,TRCBYTES);
    hdr.dt     = 1000*(dz);
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.fldr   = 1;
    hdr.trid   = TRID_DEPTH;
    hdr.ns     = n1;
    hdr.trwf   = 1;
    hdr.f1     = 0;
    hdr.f2     = 0;
    hdr.d1     = dz;
    hdr.d2     = dz;

    trrho = (float *)calloc(n1,sizeof(float));
    trvel = (float *)calloc(n1,sizeof(float));

    is = 0;
    iz = 0;
    trrho[iz]=rho[is]*1000.0;
    trvel[iz]=vel[is];
    for(iz = 1; iz < n1; iz++) {
        if (dep[is+1] > iz*dz) {
            trrho[iz]=rho[is]*1000.0;
            trvel[iz]=vel[is];
        }
        else {
            is++;
            trrho[iz]=rho[is]*1000.0;
            trvel[iz]=vel[is];
        }
        fprintf(stderr,"i=%d is=%d depth=%f iz*dz=%f rho=%f trrho=%f\n",iz, is, dep[is+1], iz*dz, rho[is], trrho[iz]);
    } 
    strcpy(filename, file_base);
    name_ext(filename, "_ro");

    fpsu = fopen(filename,"w");
    assert(fpsu != NULL);
    nwrite = fwrite( &hdr, 1, TRCBYTES, fpsu);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( &trrho[0], sizeof(float), n1, fpsu);
    assert(nwrite == n1);
    fclose(fpsu);

    strcpy(filename, file_base);
    name_ext(filename, "_cp");
    fpsu = fopen(filename,"w");
    assert(fpsu != NULL);
    nwrite = fwrite( &hdr, 1, TRCBYTES, fpsu);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( &trvel[0], sizeof(float), n1, fpsu);
    assert(nwrite == n1);
    fclose(fpsu);

    return 0;
}
