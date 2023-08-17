#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "segy.h"

/**
* writes an 2D array to a SU file
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2)
{
    size_t nwrite;
    int i;

    for (i=0; i<n2; i++) {
        nwrite = fwrite(&hdrs[i], 1, TRCBYTES, fp);
        assert(nwrite == TRCBYTES);
        nwrite = fwrite(&data[i*n1], sizeof(float), n1, fp);
        assert (nwrite == n1);
    }

    return 0;
}

#ifdef _NETCDF
#include <netcdf.h>

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int writeNetCDF(char *filename, float *data, segy *hdrs, int n1, int n2)
{
    size_t nwrite;
    int i;
    int ncid, t_dimid, x_dimid, varid;
    int dimids[2];
    size_t chunks[2];
    int shuffle, deflate, deflate_level;
    int x, y, retval;

    if (strstr(filename, ".nc") != NULL) { /* netcdf output file */
      /* Set chunking, shuffle, and deflate. */
      shuffle = NC_NOSHUFFLE;
      deflate = 1;
      deflate_level = 1;
  
      if ((retval = nc_create(filename, NC_NETCDF4, &ncid)))
      ERR(retval);
  
      /* Define the dimensions. */
      if ((retval = nc_def_dim(ncid, "t", n1, &t_dimid)))
         ERR(retval);
      if ((retval = nc_def_dim(ncid, "x", n2, &x_dimid)))
         ERR(retval);
  
      /* Set up variable data. */
      dimids[0] = t_dimid;
      dimids[1] = x_dimid;
      chunks[0] = n1;
      chunks[1] = n2/4;
  
      /* Define the variable. */
      if ((retval = nc_def_var(ncid, "data", NC_FLOAT, 2, dimids, &varid)))
         ERR(retval);
      if ((retval = nc_def_var_chunking(ncid, varid, 0, &chunks[0]))) 
         ERR(retval);
      if ((retval = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)))
         ERR(retval);
  
      /* No need to explicitly end define mode for netCDF-4 files. Write
       * the pretend data to the file. */
      if ((retval = nc_put_var_float(ncid, varid, data)))
         ERR(retval);
  
      /* Close the file. */
      if ((retval = nc_close(ncid)))
         ERR(retval);
  
      printf("*** SUCCESS writing netcdf file %s\n", filename);
    }
    else { /* writing SU file */
      FILE *fp;
      fp = fopen(filename,"w");
      assert(fp != NULL);
      for (i=0; i<n2; i++) {
        nwrite = fwrite(&hdrs[i], 1, TRCBYTES, fp);
        assert(nwrite == TRCBYTES);
        nwrite = fwrite(&data[i*n1], sizeof(float), n1, fp);
        assert (nwrite == n1);
      }
    }
    

    return 0;
}

#endif
