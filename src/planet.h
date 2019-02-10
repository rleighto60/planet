/*==============================================================================
                 ____                                                   ____
                / / /                                                  / / /
               / / /        PLANET GENERATOR                          / / /
              / / /                                                  / / /
   ____      / / /          BY: Russell Leighton          ____      / / /
   \ \ \    / / /               September 1990            \ \ \    / / /
    \ \ \  / / /                                           \ \ \  / / /
     \ \ \/ / /                                             \ \ \/ / /
      \ \/ / /                                               \ \/ / /
       \/_/_/                                                 \/_/_/

================================================================================

This program is not public domain.  It may not be distributed, sold or used in 
any commercial software without prior consent. 

Copyright © 1990,1991,1992 by Russell Leighton

=============================================================================*/

#include <stdio.h>

#include <ctype.h>
#include <malloc.h>
#include <string.h>
#define FALSE 0
#define TRUE 1
typedef unsigned char UBYTE;
typedef unsigned short UWORD;
typedef unsigned long ULONG;
typedef short WORD;
typedef char *STRPTR;

#define INCL_RXSUBCOM
#define INCL_RXSHV
#define INCL_RXFUNC
#define ULONG_TYPEDEFED
#include <rexxsaa.h>

SHVBLOCK result = {
    NULL,
    { 6, (char *)"result" },
    { 0, NULL },
    6,
    320,
    RXSHV_SYSET,
    0
};

APIRET APIENTRY planetHandler(PRXSTRING, PUSHORT, PRXSTRING);

#define RED        0x0001
#define GREEN      0x0002
#define BLUE       0x0004
#define ALPHA      0x0008
#define RMAP       0x0010
#define GMAP       0x0020
#define BMAP       0x0040

#define PLNT       0x0001
#define MOON       0x0002
#define RING       0x0004
#define PICT       0x0008

#ifndef PI
#define	PI	((double)	3.141592653589793)
#endif

char hostname[] = "PLANET";

int close_down = FALSE;
int outstanding = 0;

unsigned int mask=0,type=0,made=0;
int rc,gc,bc,rb,gb,bb;
float ra,ga,ba;
int dang,wrap,pixel;
int rad,prad,mrad,rad1,rad2,shad,olap,tran;
float trans[3][3],moon[3][3];
float light[3][3],shade[3][3];
float aspect;
int pdist,mdist,patmo;
int print=FALSE;
int pipe=FALSE;
long xres,yres,zres,vres,wres,xcen,ycen;
UBYTE *rbuf,*gbuf,*bbuf,*abuf;
UBYTE *rmap,*gmap,*bmap;

int Width,Height;

#define PGM_FORMAT ((long)'P' * 256L + (long)'2')
#define RPGM_FORMAT ((long)'P' * 256L + (long)'5')
#define PPM_FORMAT ((long)'P' * 256L + (long)'3')
#define RPPM_FORMAT ((long)'P' * 256L + (long)'6')

void *malloc();
void *realloc();
double pow(),sin(),cos(),tan(),asin(),acos(),atan2(),sqrt();
