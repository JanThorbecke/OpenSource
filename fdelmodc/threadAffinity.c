#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifdef __USE_GNU
#include <omp.h>
#include <sched.h>
#else /* for OSX */
#include <sched.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define CPU_SETSIZE	1024
#define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"


typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline void
CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
{
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int     i;
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}
#endif

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */

static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return(str);
}

void vmess(char *fmt, ...);

void threadAffinity(void)
{
  int thread;
  cpu_set_t coremask;
  char clbuf[7 * CPU_SETSIZE], hnbuf[64];
  char prefix[200];

  memset(clbuf, 0, sizeof(clbuf));
  memset(hnbuf, 0, sizeof(hnbuf));
  (void)gethostname(hnbuf, sizeof(hnbuf));

  strcpy(prefix,"Hello world from");

//  #pragma omp parallel private(thread, coremask, clbuf)
/* for use inside parallel region */
  #pragma omp critical
  {
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 1;
#endif
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
    vmess("%s thread %d, on %s. (core affinity = %s)", prefix, thread, hnbuf, clbuf);

  }
  return;
}

