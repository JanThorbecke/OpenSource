#ifndef ZFPMAR_H
#define ZFPMAR_H
#define MARBYTES        48
#define TOPBYTES        80

/* Size of the different types
long    -   8 bytes
int     -   4 bytes
char    -   1 byte
float   -   4 bytes
*/

/* TYPEDEFS */
typedef struct {	/* zfpmar - headers for the compression and decompression of data for marchenko */

    long nx; /* Number of samples in x-direction */

    long ny; /* Number of samples in y-direction */

    int sx; /* Source coordinate in x-direction */
    
    int sy; /* Source coordinate in y-direction */
    
    int sz; /* Source coordinate in z-direction */
    
    int gx; /* Receiver coordinate in x-direction */
    
    long gy; /* Receiver coordinate in y-direction */

    long compsize; /* Size of the compressed data */

} zfpmar;

typedef struct {	/* zfptop - headers for the compression and decompression of data for marchenko */
    
    long ndim; /* Number of dimension is between 1 and 4 */

    long nz; /* Number of samples in z-direction */

    long ns; /* Number of shots */

    long nt; /* Number of time samples */
    
    float dx; /* Sampling distance in x-direction */
    
    float dy; /* Sampling distance in x-direction */
    
    float dz; /* Sampling distance in x-direction */

    float scale; /* Scaling of the coordinates and the sampling */

    double tolerance; /* Set the tolerance of the zfp (de)compression */

    float fmin; /* Minimum frequency of the signal */

    float fmax; /* Maximum frequency of the singal */

    float fx; /* First location of receiver in x-direction */

    float fy; /* First location of receiver in y-direction */

    float fz; /* First location of receiver in y-direction */

} zfptop;
#endif