#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>

typedef struct _compType{ /* Receiver Type */
	int vz;  // Vertical  Particle Velocity
	int vx;  // Horizontal Particle Velocity
	int p;   // Acoustic Pressure
	int txx; // Elastic Txx
	int tzz; // Elastic Tzz
	int txz; // Elastic Txz
	int pp;  // Elastic Pressure Potential
	int ss;  // Elastic Shear    Potential
	int ud;  // ???
	int pu;  // Acoustic Up-Going    Pressure
	int pd;  // Acoustic Down-Going  Pressure
	int pl;  // Acoustic Left-Going  Pressure
	int pr;  // Acoustic Right-Going Pressure
	int pn;  // Acoustic Normal      Pressure
} compType;

typedef struct _fftPlansPar{
	// Generic Plans
	fftw_plan  fft_1d_r2r;   //One  dimensional forward real         to half-complex FFT plan
	fftw_plan ifft_1d_r2r;   //One  dimensional inverse half-complex to real         FFT plan
	fftw_plan  fft_1d_r2c;   //One  dimensional forward real         to      complex FFT plan
	fftw_plan ifft_1d_c2r;   //One  dimensional inverse      complex to real         FFT plan
	fftw_plan  fft_1d_c2c;   //One  dimensional forward      complex to      complex FFT plan
	fftw_plan ifft_1d_c2c;   //One  dimensional inverse      complex to      complex FFT plan
	fftw_plan  fft_2d_r2c;   //Two  dimensional forward real         to      complex FFT plan
	fftw_plan ifft_2d_c2r;   //Two  dimensional inverse real         to      complex FFT plan
	fftw_plan  fft_2d_c2c;   //Two  dimensional forward      complex to      complex FFT plan
	fftw_plan ifft_2d_c2c;   //Two  dimensional inverse      complex to      complex FFT plan
	fftw_plan  fft_2d_r2c_1; //Slow dimension   forward real         to      complex FFT plan
	fftw_plan ifft_2d_c2r_1; //Slow dimension   inverse real         to      complex FFT plan
	fftw_plan  fft_2d_c2c_1; //Slow dimension   forward      complex to      complex FFT plan
	fftw_plan ifft_2d_c2c_1; //Slow dimension   inverse      complex to      complex FFT plan
	fftw_plan  fft_2d_r2c_2; //Fast dimension   forward real         to      complex FFT plan
	fftw_plan ifft_2d_c2r_2; //Fast dimension   inverse real         to      complex FFT plan
	fftw_plan  fft_2d_c2c_2; //Fast dimension   forward      complex to      complex FFT plan
	fftw_plan ifft_2d_c2c_2; //Fast dimension   inverse      complex to      complex FFT plan
	///////////////////////////////////// Speacial Plans ////////////////////////////////////
	// 1D FFTw Plans
	fftw_plan  fft_1d_c2c_x;     // (x)     to (kx)    Single Fourier Transform
	fftw_plan ifft_1d_c2c_Kx;    // (kx)    to (x)     Single Fourier Transform
	fftw_plan  fft_1d_c2c_z;     // (z)     to (kz)    Single Fourier Transform
	fftw_plan ifft_1d_c2c_Kz;    // (kz)    to (z)     Single Fourier Transform
	fftw_plan  fft_1d_c2c_t;     // (t)     to (w)     Single Fourier Transform
	fftw_plan ifft_1d_c2c_W;     // (w)     to (t)     Single Fourier Transform
	// Up-Down Plane-Wave Wavefield Decomposition
	fftw_plan  fft_2d_r2c_TZ;    // (t,z)   to (w,kz)  Double Fourier Transform
	fftw_plan ifft_2d_c2r_WKz;   // (w,kz)  to (t,z)   Double Fourier Transform
	fftw_plan  fft_2d_c2c_2_WZ;  // (w,z)   to (w,kz)  Single Fourier Transform
	fftw_plan ifft_2d_c2c_2_WKz; // (w,kz)  to (w,z)   Single Fourier Transform
	//Error Next Two
	fftw_plan  fft_2d_c2c_2_TZ;  // (t,z)   to (t,kz)  Single Fourier Transform
	fftw_plan ifft_2d_c2c_2_TKz; // (t,kz)  to (t,z)   Single Fourier Transform
	// 2D Wavenumber Transform
	fftw_plan  fft_2d_r2c_ZX;    // (x,z)   to (kx,kz) Double Fourier Transform
	fftw_plan ifft_2d_c2r_KzKx;  // (kx,kz) to (x,z)   Double Fourier Transform
} fftPlansPar;

typedef struct _snapshotPar{ /* Snapshot Parameters */
	/* Filenames */
	char *file_snap; // Snapshot File Name
	char *file_beam; // Beam     File Name
	/* Snapshot Parameters */
	float tsnap1;    // First Snapshot Time
	float tsnap2;    // Last  Snapshot Time
	float dt;        // Temporal Snapshot Rate
	float dx;        // Horizontal Snapshot Sampling Rate
	float dz;        // Vertical   Snapshot Sampling Rate
	float xsnap1;    // Horizontal Coordinate Of Upper-Left  Corner Of Snapshot
	float xsnap2;    // Horizontal Coordinate Of Lower-Right Corner Of Snapshot
	float zsnap1;    // Vertical   Coordinate Of Upper-Left  Corner Of Snapshot
	float zsnap2;    // Vertical   Coordinate Of Lower-Right Corner Of Snapshot
	int ntr;
	int tracl;
	/* Snapshot Booleans */
	int forw;
	int back;
	int withbnd;
	int vxshift;
	int vzshift;
	int vxtime;
	int vztime;
	int beam;
	int decomp;
	/* Snapshot Sizes & Indices*/
	size_t nsnap;
	size_t isnap;
	size_t delay;
	size_t dtskip;
	size_t dxskip;
	size_t dzskip;
	size_t nt;
	size_t nx;
	size_t nz;
	size_t t1;
	size_t t2;
	size_t x1;
	size_t x2;
	size_t z1;
	size_t z2;
	compType type;
//	float *beam_vz;
//	float *beam_vx;
//	float *beam_p;
//	float *beam_txx;
//	float *beam_tzz;
//	float *beam_txz;
//	float *beam_pp;
//	float *beam_ss;
} snaPar;

typedef struct _modelPar{ /* Model Parameters */
	bool changedT;      //Has the time vector been changed?
	int iorder;
	int ischeme;
	int sh;
	int fldr;          // Current Field Record Number
	unsigned int dtus; // Modelling Time Step in Micro Seconds
	float dt;          // Modelling Time Step
	float dx;          // Horizontal Space Step
	float dz;          // Vertical Space Step
	float origx;       // Horizontal Model Origin (Upper Left  Corner)
	float origz;       // Vertical   Model Origin (Upper Left  Corner)
	float xmax;        // Horizontal Model End    (Lower Right Corner)
	float zmax;        // Vertical   Model End    (Lower Right Corner)
	float tmod;        // Total Modelling Time
	size_t nt;         // Number of Modelling Time Steps
	size_t nx;         // Number of Horizontal Grid Points in Model
	size_t nz;         // Number of Vertical   Grid Points in Model
	size_t nax;        // Total Number of Horizontal Grid Points
	size_t naz;        // Total Number of Vertical   Grid Points
	size_t sizem;      // Total Modelling Domain Size
	/* Vx: rox */
	size_t ioXxb;      // Horizontal Boundary Layer Start
	size_t ieXxb;      // Horizontal Boundary Layer End
	size_t ioXzb;      // Vertical   Boundary Layer Start
	size_t ieXzb;      // Vertical   Boundary Layer End
	size_t ioXx;       // Horizontal Model    Layer Start
	size_t ieXx;       // Horizontal Model    Layer End
	size_t ioXz;       // Vertical   Model    Layer Start
	size_t ieXz;       // Vertical   Model    Layer End
	/* Vz: roz */
	size_t ioZxb;      // Horizontal Boundary Layer Start
	size_t ieZxb;      // Horizontal Boundary Layer End
	size_t ioZzb;      // Vertical   Boundary Layer Start
	size_t ieZzb;      // Vertical   Boundary Layer End
	size_t ioZx;       // Horizontal Model    Layer Start
	size_t ieZx;       // Horizontal Model    Layer End
	size_t ioZz;       // Vertical   Model    Layer Start
	size_t ieZz;       // Vertical   Model    Layer End
	/* P, Txx, Tzz: lam, l2m */
	size_t ioPxb;      // Horizontal Boundary Layer Start
	size_t iePxb;      // Horizontal Boundary Layer End
	size_t ioPzb;      // Vertical   Boundary Layer Start
	size_t iePzb;      // Vertical   Boundary Layer End
	size_t ioPx;       // Horizontal Model    Layer Start
	size_t iePx;       // Horizontal Model    Layer End
	size_t ioPz;       // Vertical   Model    Layer Start
	size_t iePz;       // Vertical   Model    Layer End
	/* Txz: muu */
	size_t ioTxb;      // Horizontal Boundary Layer Start
	size_t ieTxb;      // Horizontal Boundary Layer End
	size_t ioTzb;      // Vertical   Boundary Layer Start
	size_t ieTzb;      // Vertical   Boundary Layer End
	size_t ioTx;       // Horizontal Model    Layer Start
	size_t ieTx;       // Horizontal Model    Layer End
	size_t ioTz;       // Vertical   Model    Layer Start
	size_t ieTz;       // Vertical   Model    Layer End
	/* Model Statistics */
	float cp_max;
	float cp_min;
	float cs_max;
	float cs_min;
	float rho_max;
	float rho_min;
	/* Filenames */
	char *file_cp;  //Acoustic Velocity
	char *file_cs;  //Shear Wave Velocity
	char *file_den; //Density
	char *file_qp;  //
	char *file_qs;  //
	char *file_imp; //Impedance
	char *file_dd;  //Decomposition Direction
	/* Models */
	float *cp;  /* P-Wave Velocity Model */
	float *cs;  /* S-Wave Velocity Model */
	float *rho; /* Density         Model */
	float *qp;
	float *qs;
	float *rox;
	float *roz;
	float *l2m;
	float *tss;
	float *tes;
	float *tep;
	float *lam;
	float *mul;
	float *r;
	float *p;
	float *q;
	float *imp;   /*                                  Acoustic Impedance          */
	float *ngxv;  /* Horizontal            Normalized Acoustic Impedance Gradient */
	float *ngzv;  /* Vertical              Normalized Acoustic Impedance Gradient */
//	float *nogxv; /* Horizontal Orthogonal Normalized Acoustic Impedance Gradient */
//	float *nogzv; /* Vertical   Orthogonal Normalized Acoustic Impedance Gradient */
	/*???*/
	float z0;
	float x0;
	float Qp;
	float Qs;
	float fw;
	int mavgi; //Impedance Model Moving Average Filter
} modPar;

typedef struct _wavfieldPar{ /* Wavefield Paramters */
	float *tzz;   // Pressure (acoustic)                      Wavefield
	float *txz;  //
	float *txx;  //
	float *r;    //
	float *p;    //
	float *q;    //
	float *vx;   // Horizontal  Particle            Velocity Wavefield
	float *vz;   // Vertical    Particle            Velocity Wavefield
	float *dtzz; // Tzz/P       Temporal                     Gradient
	float *dvx;  // Horizontal  Particle            Velocity Gradient
	float *dvz;  // Vertical    Particle            Velocity Gradient
	float *pu;   // Up-Going    Pressure                     Wavefield
	float *pd;   // Down-Going  Pressure                     Wavefield
	float *pl;   // Left-Going  Pressure                     Wavefield
	float *pr;   // Right-Going Pressure                     Wavefield
	float *pn;   // Normal      Pressure                     Wavefield
//	float *vxu;  // Up-Going    Horizontal Particle Velocity Wavefield
//	float *vxd;  // Down-Going  Horizontal Particle Velocity Wavefield
//	float *vxl;  // Left-Going  Horizontal Particle Velocity Wavefield
//	float *vxr;  // Right-Going Horizontal Particle Velocity Wavefield
//	float *vzu;  // Up-Going    Vertical   Particle Velocity Wavefield
//	float *vzd;  // Down-Going  Vertical   Particle Velocity Wavefield
//	float *vzl;  // Left-Going  Vertical   Particle Velocity Wavefield
//	float *vzr;  // Right-Going Vertical   Particle Velocity Wavefield
} wavPar;

typedef struct _sourcePar{ /* Source Array Parameters */
	char *file_src; //Source Filename
	int *typ;       // Source Type
	int *orient;    // Source Orientation
	size_t *ind;    // Index      Source Grid Point
	size_t *xi;     // Horizontal Source Grid Point
	size_t *zi;     // Vertical   Source Grid Point
	float *x;       // Horizontal Source Location
	float *z;       // Source Depth
	float *wav;     // Wavefield (P, Vx or Vz)
	off_t loc;      // Current position of file pointer
	bool eof;       // End Of File (EOF) Reached?
	int tracl;      // Current Trace Number In File
	size_t nsrc;    // Number of Sources in Current Field Record
	size_t ntrc;    // Number of Processed Traces
	size_t nt;      // Number of Time Samples per Trace
} srcPar;

typedef struct _receiverPar{ /* Receiver Parameters */
	size_t nrcv;    // Number of Receivers
	char *file_rcv; // Receiver Filename
	int *typ;       // Receiver Type
	int *orient;    // Receiver Orientation
	size_t *loc;    // Receiver Index      Location in Model
	size_t *xi;     // Horizontal Receiver Location in Model
	size_t *zi;     // Vertical   Receiver Location in Model
	float *wav;     // Stored Wavefield
} rcvPar;

typedef struct _shotPar{ /* Shot Parameters */
	int n;
	int *z;
	int *x;
} shotPar;

typedef struct _boundPar{ /* Boundary Parameters */
	int top;
	int bot;
	int lef;
	int rig;
	int cfree;
	int pml;
	size_t ntap;
	size_t ntapo;
	size_t npml;
	float R;
	float m;
	float tapfact;
	size_t *surface;
	float *tapz;
	float *tapx;
	float *tapxz;
	float *pml_Vx;
	float *pml_nzVx;
	float *pml_nxVz;
	float *pml_nzVz;
	float *pml_nxP;
	float *pml_nzP;
} bndPar;

typedef struct _decompPar{ /* Decomposition Parameters */
	// Decomposition Booleans
	int decomp;  //Boolean - Decomposition If True
	int direct;  //Integer - 0: Nothing 1: Down-Going 2: Up-Going 3: Left-Going 4: Right-Going
	int pu;      //Boolean - Up-Going Pressure
	int pd;      //Boolean - Down-Going Pressure
	int pr;      //Boolean - Right-Going Pressure
	int pl;      //Boolean - Left-Going Pressure
	int pn;      //Boolean - Normal-Up-Going Pressure
	int vxu;     //Boolean - Up-Going Horizontal Particle Velocity
	int vxd;     //Boolean - Down-Going Horizontal Particle Velocity
	int vxr;     //Boolean - Right-Going Horizontal Particle Velocity
	int vxl;     //Boolean - Left-Going Horizontal Particle Velocity
	int vxn;     //Boolean - Normal-Up-Going Horizontal Particle Velocity
	int vzu;     //Boolean - Up-Going Vertical Particle Velocity
	int vzd;     //Boolean - Down-Going Vertical Particle Velocity
	int vzr;     //Boolean - Right-Going Vertical Particle Velocity
	int vzl;     //Boolean - Left-Going Vertical Particle Velocity
	int vzn;     //Boolean - Normal-Up-Going Vertical Particle Velocity
	int wavFilt; //Boolean - Wavenumber Filter
	// Regularization & Decomposition Operator
	int med;     //Median Order 3 or 5;
	int mavgn;   //Integer - Moving Average Preferential Direction Filter Order 0,1,3,5,7,9
	int mavgo;   //Integer - Moving Average Orthogonal   Direction Filter Order 0,1,3,5,7,9
	int mavga;   //Integer - Moving Average Decomp. Dir. Angle     Filter Order 0,1,3,5,7,9
	int writeDD; //Boolean - Write Out Decomposition Direction?
	float reg;   //Regularization Parameter
	float px;    //Horizontal Magnitude of Preferential Direction Vector
	float pz;    //Vertical   Magnitude of Preferential Direction Vector
	float kl;    //Star of Wavenumber Filter
	float kh;    //End   of Wavenumber Filter
	float *op;   //Square Root Decomposition Operator
	float *tmp;  //Temporary Storage Array
} decompPar;

typedef struct _migPar{
	int mode;             //Migration Mode 1: Conventional 2: Poynting 3: Decomposition
	int orient;           //Migration Orientation
	int backscatter;      //Image Back-Scattered Wavefields
	int compress;         //Compress Forward Modelled Snapshots?
	// Migration Image Size & Sampling Rates
	float xmig1;          //Horizontal Coordinate of Upper-Left  Corner of Migration Image
	float xmig2;          //Horizontal Coordinate of Lower-Right Corner of Migration Image
	float zmig1;          //Vertical   Coordinate of Upper-Left  Corner of Migration Image
	float zmig2;          //Vertical   Coordinate of Lower-Right Corner of Migration Image
	float dt;             //Migration Image Temporal   Sampling Rate
	float dx;             //Migration Image Horizontal Sampling Rate
	float dz;             //Migration Image Vertical   Sampling Rate
	size_t nt;            //Number of Migration Image Time Samples
	size_t nx;            //Number of Horizontal Migration Image Grid Points
	size_t nz;            //Number of Vertical   Migration Image Grid Points
	size_t sizem;         //Total Number of Grid Points in Migration Image
	size_t x1;            //Horizontal Model Grid Point of Upper-Left  Corner of Migration Image
	size_t x2;            //Horizontal Model Grid Point of Lower-Right Corner of Migration Image
	size_t z1;            //Vertical   Model Grid Point of Upper-Left  Corner of Migration Image
	size_t z2;            //Vertical   Model Grid Point of Lower-Right Corner of Migration Image
	size_t it;            //Current Migration Time Step (Multiple of mod.it)
	size_t skipdt;        //Number of Modelling Time Samples Between Imaging Condition Application
	size_t skipdx;        //Number of Horizontal Model Grid Points Between Horizontal Migration Image Grid Points
	size_t skipdz;        //Number of Vertical   Model Grid Points Between Vertical   Migration Image Grid Points
	// Output Variable
	size_t ntr;           //Cumulative Number Of Traces in "file_image"
	// Pointers
	char *file_mig;       //Filename for Storage of Migrated Image(s)
	float *image;         //Migration Image for Current fldr
	float *mig;           //Migration Image
	wavPar *wav;          //Wavefield Snapshots for Current fldr
//	decompWav *decompWav; //Decomposed Forward Propagated Snapshots
	// Plane Wave Decomposition
	int tap;     //Size of Plane-Wave Taper
	// Booleans For Decomposition Direction Combinations For Imaging
	         // (S) - (R)          S:Source Wavefield R:Receiver Wavefield
	bool pp; // (+) - (+) Imaging (Transmission)
	bool pm; // (+) - (-) Imaging ( Reflection )
	bool mp; // (-) - (+) Imaging ( Reflection )
	bool mm; // (+) - (+) Imaging (Transmission)
// NOTE: We always image the pressure, we should add the option to image the particle velocity
} migPar;

typedef struct _recPar{
	bool rec;       // Record Wavefield
	int  left;      // Left   Boundary?
	int  top;       // Top    Boundary?
	int  bottom;    // Bottom Boundary?
	int  right;     // Right  Boundary?
	int  p;         // Record Pressure?
	int  txx;       // Record Txx?
	int  txz;       // Record Txz?
	int  tzz;       // Record Tzz?
	int  vx;        // Record Vx?
	int  vz;        // Record Vz?
	int  write;     // Write Out Recorded Wavefield?
	char *file_loc; //Filename of Receiver Locations File
} recPar;


#if __STDC_VERSION__ >= 199901L
	/* "restrict" is a keyword */
#else
#define restrict 
#endif

#ifndef WISDOMDIR
#define WISDOMDIR "/tmp/fftw/"
#endif
