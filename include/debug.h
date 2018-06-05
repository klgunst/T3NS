#ifndef DEBUG_H
# define DEBUG_H
# ifndef DEBUG
#  define NDEBUG
# endif
#include <assert.h> /* asserts can be disabled by defining a macro NDEBUG, removes automatically 
                       the asserts when compiling non-debug code */
#endif
