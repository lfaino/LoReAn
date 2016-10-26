static char rcsid[] = "$Id: bigendian.c 168395 2015-06-26 17:13:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bigendian.h"
#include <unistd.h>		/* For read() */


/* Same as Littleendian_write_char */
void
Bigendian_write_char (unsigned char value, int fd) {
  unsigned char buf[1];

  buf[0] = value;
  write(fd,buf,1);

  return;
}

/************************************************************************
 *   Int
 ************************************************************************/

int
Bigendian_convert_int (int littleendian) {
  int bigendian;

  bigendian = littleendian & 0xff; /* 0 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */

  return bigendian;
}


size_t
Bigendian_fwrite_int (int value, FILE *fp) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}


size_t
Bigendian_fwrite_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    value = array[i];
    buf[3] = value & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}


size_t
Bigendian_fread_int (int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
#if 0
    *value = buf[0];
    *value <<= 8;
    *value |= buf[1];
    *value <<= 8;
    *value |= buf[2];
    *value <<= 8;
    *value |= buf[3];
#else
    *value = ((int) buf[0] << 24) | ((int) buf[1] << 16) | ((int) buf[2] << 8) | (int) buf[3];
#endif
    return 1;
  }
}


size_t
Bigendian_fread_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  /* int value; */
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
#if 0
      value = buf[0];
      value <<= 8;
      value |= buf[1];
      value <<= 8;
      value |= buf[2];
      value <<= 8;
      value |= buf[3];
      array[i] = value;
#else
      array[i] = ((int) buf[0] << 24) | ((int) buf[1] << 16) | ((int) buf[2] << 8) | (int) buf[3];
#endif
    }
  }
  return n;
}


/************************************************************************
 *   Unsigned int
 ************************************************************************/

unsigned int
Bigendian_convert_uint (unsigned int littleendian) {
  unsigned int bigendian;

  bigendian = littleendian & 0xff; /* 0 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */

  return bigendian;
}


size_t
Bigendian_fwrite_uint (unsigned int value, FILE *fp) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}


void
Bigendian_write_uint (unsigned int value, int fd) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  write(fd,buf,4);
  return;
}


size_t
Bigendian_fwrite_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  unsigned int value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[3] = value & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}


size_t
Bigendian_fread_uint (unsigned int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
#if 0
    *value = buf[0];
    *value <<= 8;
    *value |= buf[1];
    *value <<= 8;
    *value |= buf[2];
    *value <<= 8;
    *value |= buf[3];
#else
    *value = ((unsigned int) buf[0] << 24) | ((unsigned int) buf[1] << 16) | ((unsigned int) buf[2] << 8) | (unsigned int) buf[3];
#endif
    return 1;
  }
}


size_t
Bigendian_fread_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  /* unsigned int value; */
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
#if 0
      value = buf[0];
      value <<= 8;
      value |= buf[1];
      value <<= 8;
      value |= buf[2];
      value <<= 8;
      value |= buf[3];
      array[i] = value;
#else
      array[i] = ((unsigned int) buf[0] << 24) | ((unsigned int) buf[1] << 16) | ((unsigned int) buf[2] << 8) | (unsigned int) buf[3];
#endif
    }
  }
  return n;
}


unsigned int
Bigendian_fileio_read_uint (int fd) {
  unsigned int value = 0U;
  unsigned char buf[4];

  read(fd,buf,4);
#if 0
  value = buf[0];
  value <<= 8;
  value |= buf[1];
  value <<= 8;
  value |= buf[2];
  value <<= 8;
  value |= buf[3];
#else
  value = ((unsigned int) buf[0] << 24) | ((unsigned int) buf[1] << 16) | ((unsigned int) buf[2] << 8) | (unsigned int) buf[3];
#endif
  return value;
}


/************************************************************************
 *   Long unsigned int
 ************************************************************************/

#ifdef HAVE_64_BIT

UINT8
Bigendian_convert_uint8 (UINT8 littleendian) {
  UINT8 bigendian;
  unsigned char byte1, byte2, byte3, byte4, byte5, byte6, byte7;

  bigendian = littleendian & 0xff;
  byte1 = (littleendian >>= 8);
  byte2 = (littleendian >>= 8);
  byte3 = (littleendian >>= 8);
  byte4 = (littleendian >>= 8);
  byte5 = (littleendian >>= 8);
  byte6 = (littleendian >>= 8);
  byte7 = (littleendian >>= 8);

  /* bigendian = byte0; */
  bigendian <<= 8;
  bigendian |= (byte1 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte2 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte3 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte4 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte5 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte6 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte7 & 0xff);

  return bigendian;
}


void
Bigendian_write_uint8 (UINT8 value, int fd) {
  unsigned char buf[8];

  buf[7] = value & 0xff;
  buf[6] = (value >>= 8) & 0xff;
  buf[5] = (value >>= 8) & 0xff;
  buf[4] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  write(fd,buf,8);
  return;
}


size_t
Bigendian_fwrite_uint8 (UINT8 value, FILE *fp) {
  unsigned char buf[8];

  buf[7] = value & 0xff;
  buf[6] = (value >>= 8) & 0xff;
  buf[5] = (value >>= 8) & 0xff;
  buf[4] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}


size_t
Bigendian_fwrite_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[7] = value & 0xff;
    buf[6] = (value >>= 8) & 0xff;
    buf[5] = (value >>= 8) & 0xff;
    buf[4] = (value >>= 8) & 0xff;
    buf[3] = (value >>= 8) & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}


size_t
Bigendian_fread_uint8 (UINT8 *value, FILE *fp) {
  unsigned char buf[8];

  if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (UINT8) buf[0];
    *value <<= 8;
    *value |= (UINT8) buf[1];
    *value <<= 8;
    *value |= (UINT8) buf[2];
    *value <<= 8;
    *value |= (UINT8) buf[3];
    *value <<= 8;
    *value |= (UINT8) buf[4];
    *value <<= 8;
    *value |= (UINT8) buf[5];
    *value <<= 8;
    *value |= (UINT8) buf[6];
    *value <<= 8;
    *value |= (UINT8) buf[7];
    return 1;
  }
}


size_t
Bigendian_fread_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (UINT8) buf[0];
      value <<= 8;
      value |= (UINT8) buf[1];
      value <<= 8;
      value |= (UINT8) buf[2];
      value <<= 8;
      value |= (UINT8) buf[3];
      value <<= 8;
      value |= (UINT8) buf[4];
      value <<= 8;
      value |= (UINT8) buf[5];
      value <<= 8;
      value |= (UINT8) buf[6];
      value <<= 8;
      value |= (UINT8) buf[7];
      array[i] = value;
    }
  }
  return n;
}


UINT8
Bigendian_fileio_read_uint8 (int fd) {
  UINT8 value = 0LU;
  unsigned char buf[8];

  read(fd,buf,8);
  value = (UINT8) buf[0];
  value <<= 8;
  value |= (UINT8) buf[1];
  value <<= 8;
  value |= (UINT8) buf[2];
  value <<= 8;
  value |= (UINT8) buf[3];
  value <<= 8;
  value |= (UINT8) buf[4];
  value <<= 8;
  value |= (UINT8) buf[5];
  value <<= 8;
  value |= (UINT8) buf[6];
  value <<= 8;
  value |= (UINT8) buf[7];
  return value;
}

#endif /* HAVE_64_BIT */


/************************************************************************
 *   Double
 ************************************************************************/

size_t
Bigendian_fwrite_double (double value, FILE *fp) {
  unsigned char buf[8], *ptr = (unsigned char *) &value;
  size_t i, j;

  /* buf = (unsigned char *) MALLOC(sizeof(double) * sizeof(unsigned char)); */

  i = 0;
  j = sizeof(double);
  while (i < sizeof(double)) {
    buf[i++] = ptr[--j];
  }

  if (fwrite(buf,sizeof(unsigned char),sizeof(double),fp) == 0) {
    /* Should set error indicator for stream and set errno */
    /* FREE(buf); */
    return 0;
  } else {
    /* FREE(buf); */
    return sizeof(double)/4;
  }
}


double
Bigendian_convert_double (double value) {
  unsigned char *ptr = (unsigned char *) &value, temp;
  size_t i, j;

  i = 0;
  j = sizeof(double);
  while (i < sizeof(double)) {
    /* swap */
    temp = ptr[--j];
    ptr[j] = ptr[i];
    ptr[i++] = temp;
  }

  return value;
}

