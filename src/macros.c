#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "macros.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_status(void);

/* ========================================================================== */

void* safe_malloc_helper(long long s, size_t t, const char *typ, const char *file, int line, const char *func){

  void *pn = malloc(s * t);
  if (pn == NULL || s < 0){
    const char *filname = strrchr(file, '/') + 1;
    fprintf(stderr, "%s:%d @%s :: Failed to allocate %s array of size %lld (%llu bytes)!\n"
            "Maximal size of size_t : %ld\n",
        filname == NULL ? file : filname, line, func, typ, s, s*t, SIZE_MAX);
    print_status();
    exit(EXIT_FAILURE);
  }
  return pn;
}
 
void* safe_calloc_helper(long long s, size_t t, const char *typ, const char *file, int line, const char *func){

  void *pn = calloc(s, t);
  if (pn == NULL || s < 0){
    const char *filname = strrchr(file, '/') + 1;
    fprintf(stderr, "%s:%d @%s :: Failed to reallocate %s array of size %lld (%llu bytes)!\n"
            "Maximal size of size_t : %ld\n",
        filname == NULL ? file : filname, line, func, typ, s, s*t, SIZE_MAX);
    print_status();
    exit(EXIT_FAILURE);
  }
  return pn;
}

void safe_free_helper(void *ptr, const char *ptrname, const char *file, int line, const char *func){

  if (ptr == NULL){
    const char *filname = strrchr(file, '/') + 1;
    fprintf(stderr, "%s:%d @%s :: Pointer %s is a NULL-pointer, can not be freed!\n",
        filname == NULL ? file : filname, line, func, ptrname);
    exit(EXIT_FAILURE);
  }
  safe_free(ptr);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void print_status(void){
  /* /proc/[pid]/status contains all human readible status, if i want to get the memory later 
   * on and am only interested in that (e.g. to write to a file every x seconds), 
   * it is better to read from /proc/[pid]/statm (easier to parse)
   */
  FILE* status;
  char line[256];
  if ((status = fopen("/proc/self/status", "r")) == NULL){
    fprintf(stderr, "Error in openeing /proc/self/status\n");
    return;
  }
  while (fgets(line, sizeof line, status)){
    fprintf(stderr, "%s", line);
  }
  fclose(status);
}
