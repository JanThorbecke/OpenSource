#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef WAVEPAR
#define WAVEPAR
typedef struct WaveParameters {
    int nt, shift, inv, scfft, cm, cn, wav;
    float dt, fp, fmin, flef, frig, fmax, t0, db, scale, eps;
    char w[10], *file_wav;
} WavePar;
#endif

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

#ifndef FQWV
#define FQWV
void freqwave(float *wave, int nt, float dt, float fp, float fmin, float flef, float frig, float fmax, float t0, float db, int shift, int cm, int cn, char *w, float scale, int scfft, int inverse, float eps, int verbose);
#endif
