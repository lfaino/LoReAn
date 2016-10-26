#ifndef FOPEN_INCLUDED
#define FOPEN_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>		/* For USE_FOPEN_BINARY, USE_FOPEN_TEXT */
#endif


#if USE_FOPEN_BINARY
#define FOPEN_READ_BINARY(filename) fopen(filename,"rb")
#define FOPEN_WRITE_BINARY(filename) fopen(filename,"wb")
#define FOPEN_RW_BINARY(filename) fopen(filename,"w+b")
#else
#define FOPEN_READ_BINARY(filename) fopen(filename,"r")
#define FOPEN_WRITE_BINARY(filename) fopen(filename,"w")
#define FOPEN_RW_BINARY(filename) fopen(filename,"w+")
#endif

#if USE_FOPEN_TEXT
#define FOPEN_READ_TEXT(filename) fopen(filename,"rt")
#define FOPEN_WRITE_TEXT(filename) fopen(filename,"wt")
#define FOPEN_RW_TEXT(filename) fopen(filename,"w+t")
#else
#define FOPEN_READ_TEXT(filename) fopen(filename,"r")
#define FOPEN_WRITE_TEXT(filename) fopen(filename,"w")
#define FOPEN_RW_TEXT(filename) fopen(filename,"w+")
#endif


#endif

