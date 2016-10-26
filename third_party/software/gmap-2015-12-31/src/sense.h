#ifndef SENSE_INCLUDED
#define SENSE_INCLUDED

#define SENSE_NULL 0x0
#define SENSE_ANTI 0x1
#define SENSE_FORWARD 0x2

#define SENSE_CONSISTENT_P(x,y) ((x | y) != 0x3)
#define SENSE_INCONSISTENT_P(x,y) ((x | y) == 0x3)

#define SENSE_CONSISTENT_FOR_INVERSION_P(x,y) ((x & y) == 0x0)
#define SENSE_INCONSISTENT_FOR_INVERSION_P(x,y) ((x & y) != 0x0)


#endif

