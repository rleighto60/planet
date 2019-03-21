/*==============================================================================
                 ____                                                   ____
                / / /                                                  / / /
               / / /        PLANET GENERATOR v3.4                     / / /
              / / /                                                  / / /
   ____      / / /          BY: Russell Leighton          ____      / / /
   \ \ \    / / /               August 1991               \ \ \    / / /
    \ \ \  / / /                                           \ \ \  / / /
     \ \ \/ / /                                             \ \ \/ / /
      \ \/ / /                                               \ \/ / /
       \/_/_/                                                 \/_/_/

================================================================================

This program is not public domain.  It may not be distributed, sold or used in 
any commercial software without prior consent. 

Copyright © 1990,1991,1992 by Russell Leighton

==============================================================================*/

#include "planet.h"

extern int parse(char *input, char *args[], int narg);

/* Equate two matrices M = N */

void mset(float (*n)[3], float (*m)[3])
{
    int i,j;

    for(j = 0; j < 3; j++) {
        for(i = 0; i < 3; i++) m[i][j] = n[i][j];
    }
}

/* Multiply two matrices together M = N*M */

void mmul(float (*n)[3], float (*m)[3])
{
    int i,j,k;
    float p[3][3];

    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            p[i][j] = 0.0;
            for(k = 0; k < 3; k++) p[i][j] += n[i][k]*m[k][j];
        }
    }
    mset(p,m);
}

/* compute transformation and concatenate with previous transformation(s) */

void rotate(int axis, int ang, float (*m)[3])
{
   int i,j,k;
   float n[3][3];

   if(!ang) return;

/* setup tranformation matrix based on selected axis */

   switch(axis) {
      case 'x': i = 0; j = 1; k = 2; break;
      case 'y': i = 1; j = 2; k = 0; break;
      case 'z': i = 2; j = 0; k = 1; break;
      default : i = 0; j = 1; k = 2; break;
   }

   n[i][i] = 1.0;
   n[j][j] = n[k][k] = (float)cos((double)ang*PI/180.0);
   n[i][j] = n[i][k] = n[j][i] = n[k][i] = 0.0;
   n[k][j] = (float)sin((double)ang*PI/180.0);
   n[j][k] = -n[k][j];

   mmul(n,m);

} /* end rotate */

/* Return identity matrix M = I */

void identity(float (*m)[3])
{
   int i,j;

   for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
         if(i == j) m[i][j] = 1;
         else m[i][j] = 0;
      }
   }
}

/* get color from map at specified px,py,pz */

void map(long vx, long vy, long vz, long va, float ca, float f)
{
   long x,y;
   long pos;
   long v;
   int r1,g1,b1,r2,g2,b2,c;
   float lon,lat,cf;

/* compute longitude and latitude */

   lon = (float)atan2((double)vx,(double)(2.0-f)*vz);
   if(wrap == 1) {
      lon *= 2.0;
      if(lon < 0.0) lon += PI;
      else lon -= PI;
   }
   v = (float)sqrt((double)(vx*vx + vz*vz));
   lat = (float)atan2((double)f*vy,(double)v);

   x = (.5 + lon/(2.0*PI))*(Width-olap-1);
   y = (.5 - lat/PI)*(Height-1);

   pos = y*Width+x;
   r1 = (int)rmap[pos];
   g1 = (int)gmap[pos];
   b1 = (int)bmap[pos];
   if((x < olap) && (olap != 0)) {
      pos += Width-olap;
      r2 = (int)rmap[pos];
      g2 = (int)gmap[pos];
      b2 = (int)bmap[pos];
      r1 = (x*r1 + (olap-x)*r2)/olap;
      g1 = (x*g1 + (olap-x)*g2)/olap;
      b1 = (x*b1 + (olap-x)*b2)/olap;
   }

   if((r1 >= rb) || (g1 >= gb) || (b1 >= bb)) {
      if(ca > ra) r2 = ca*r1;
      else r2 = ra*r1;
      if(ca > ga) g2 = ca*g1;
      else g2 = ga*g1;
      if(ca > ba) b2 = ca*b1;
      else b2 = ba*b1;
      if (va == 0) { rc = r2; gc = g2; bc = b2; }
      else {
         if (va > 0) {
            c = (r1 > g1) ? r1 : g1;
            c = (c > b1) ? c : b1;
            cf = (float)(c*(255-va))/65025.0;
         }
         else
            cf = 1.0 + (float)va/255.0;
         rc = cf*r2 + (1.0-cf)*rc;
         gc = cf*g2 + (1.0-cf)*gc;
         bc = cf*b2 + (1.0-cf)*bc;
      }
      pixel = TRUE;
   }
} /* end map */

void ringmap(float r, float ang, float ca)
{
   long x,y;
   long pos;
   int r1,g1,b1,r2,g2,b2,c;
   float cf;

   x = (ang/PI + 1.0)*(Width-1)/2.0;
   y = r*(Height-1);

   pos = y*Width+x;
   r1 = (int)rmap[pos];
   g1 = (int)gmap[pos];
   b1 = (int)bmap[pos];
   if(x < olap) {
      pos += Width-olap;
      r2 = (int)rmap[pos];
      g2 = (int)gmap[pos];
      b2 = (int)bmap[pos];
      r1 = (x*r1 + (olap-x)*r2)/olap;
      g1 = (x*g1 + (olap-x)*g2)/olap;
      b1 = (x*b1 + (olap-x)*b2)/olap;
   }

   if((r1 >= rb) || (g1 >= gb) || (b1 >= bb)) {
      if(ca > ra) r2 = ca*r1;
      else r2 = ra*r1;
      if(ca > ga) g2 = ca*g1;
      else g2 = ga*g1;
      if(ca > ba) b2 = ca*b1;
      else b2 = ba*b1;
      if (tran == 0) { rc = r2; gc = g2; bc = b2; }
      else {
         if (tran > 0) {
            c = (r1 > g1) ? r1 : g1;
            c = (c > b1) ? c : b1;
            cf = (float)(c*(255-tran))/65025.0;
         }
         else
            cf = 1.0 + (float)tran/255.0;
         rc = cf*r2 + (1.0-cf)*rc;
         gc = cf*g2 + (1.0-cf)*gc;
         bc = cf*b2 + (1.0-cf)*bc;
      }
      pixel = TRUE;
   }
} /* end ringmap */

int open_things()
{
   long buflen;

   vres = xres;
   if (dang) vres = 2 * xres;
   buflen = vres * yres;

   if(!(mask & RED)) {
      if ((rbuf = calloc((size_t)buflen,1L)) == NULL) {
         fprintf(stderr,"main - insufficient memory!!!\n");
         return(0);
      }
      mask |= RED;
   }
   if(!(mask & GREEN)) {
      if ((gbuf = calloc((size_t)buflen,1L)) == NULL) {
         fprintf(stderr,"main - insufficient memory!!!\n");
         return(0);
      }
      mask |= GREEN;
   }
   if(!(mask & BLUE)) {
      if ((bbuf = calloc((size_t)buflen,1L)) == NULL) {
         fprintf(stderr,"main - insufficient memory!!!\n");
         return(0);
      }
      mask |= BLUE;
   }
   buflen = vres/8L * yres;
   if(!(mask & ALPHA)) {
      if ((abuf = calloc((size_t)buflen,1L)) == NULL) {
         fprintf(stderr,"main - insufficient memory!!!\n");
         return(0);
      }
      mask |= ALPHA;
   }
   return(1);
}

void close_things()
{
   if(mask & RED) {
      free(rbuf);
      mask ^= RED;
   }
   if(mask & GREEN) {
      free(gbuf);
      mask ^= GREEN;
   }
   if(mask & BLUE) {
      free(bbuf);
      mask ^= BLUE;
   }
   if(mask & ALPHA) {
      free(abuf);
      mask ^= ALPHA;
   }
}

int generate(int lflag, int depth)
{
    int d,r=0,r1=0,r2=0;
    long pcx=0,pcy=0;
    long x,y,bpos,mpos;
    long vx,vy,vz,va,px,py,pz,sx,sy,sz,s,v,pa,pz2,rr=0,rr1=0,rr2=0,pr2;
    long mcx,mcy,mcz=0;
    long minx,maxx,miny,maxy,xoff;
    float f,calp,fdist,alp,pr,dr;
    float matr[3][3];

/* determine minimum possible window to conserve computation */

    if(type & PICT) {
         minx = 0L;
         miny = 0L;
         maxx = xres;
         maxy = yres;
    }
    else {
        if(type & MOON) {
            mset(moon,matr);
	    if(dang) {
	      if(lflag) rotate('y',-(dang>>1),matr);
	      else rotate('y',(dang>>1),matr);
	    }
            s = (mdist*xres)/wres;
            mcx = s*matr[0][2]*aspect;
            mcy = s*matr[1][2];
            mcz = s*matr[2][2];
        }
        else mcx = mcy = mcz = 0;
        d = (pdist*xres)/wres - mcz;
        if(d <= 0) return(0);
        fdist = (float)zres/(float)d;
        r = (rad*xres*fdist)/wres;
	pcx = (xcen*xres)/wres+xres/2;
	if(dang) {
	  if(lflag) pcx -= depth>>1;
	  else pcx += depth>>1;
	}
        pcy = (ycen*xres)/wres+yres/2;
        pcx += mcx*(float)wres/(float)pdist;
        pcy -= mcy*(float)wres/(float)pdist;
        if(type & RING) {
           r1 = (rad1*xres*fdist)/wres;
           r2 = (rad2*xres*fdist)/wres;
           rr1 = (long)r1*r1;
           rr = r*r;
        }
        else r2 = r;
        rr2 = (long)r2*r2;

        minx = (long)pcx - (long)(r2*aspect);
        minx = (minx > 0L) ? minx : 0L;
        miny = (long)pcy - (long)r2;
        miny = (miny > 0L) ? miny : 0L;
        maxx = (long)pcx + (long)(r2*aspect);
        maxx = (maxx < xres) ? maxx : xres;
        maxy = (long)pcy + (long)r2;
        maxy = (maxy < yres) ? maxy : yres;
        identity(matr);
	if(dang) {
	  if(lflag) rotate('y',(dang>>1),matr);
	  else rotate('y',-(dang>>1),matr);
	}
        mmul(trans,matr);
    }
    pa = (patmo*xres)/wres;

    if (lflag) xoff = xres;
    else xoff = 0;

/* loop on destination display coordinates */

   for (y = miny; y < maxy; y++) {
/*      if(lflag == (y%2)) continue; */
      if(type & PICT) {
         py = (y*(long)Height)/yres;
      }
      else py = (long)pcy - y;

      for (x = minx; x < maxx; ++x) {

         pixel = FALSE;
         bpos = vres*y+(x+xoff);
         rc = rbuf[bpos];
         gc = gbuf[bpos];
         bc = bbuf[bpos];
         dr = pa;

         switch(type) {
            case PICT:
            px = x;
	    if(dang) {
	      if(lflag) px += depth>>1;
	      else px -= depth>>1;
	    }
	    px *= (long)Width/xres;
            if((px < 0) || (px >= Width)) { rc = gc = bc = 0; break; }
            mpos = py*Width+px;
            rc = rmap[mpos];
            gc = gmap[mpos];
            bc = bmap[mpos];
            if((rc >= rb) || (gc >= gb) || (bc >= bb)) pixel = TRUE;
            break;

            case RING:
            px = (x - (long)pcx)/aspect;
            if(matr[1][2] == 0.0) break;
            pz = -(px*matr[1][0] + py*matr[1][1])/matr[1][2];
            if(pz < 0) {
	      if((abuf[bpos>>3] >> (7-((x+xoff) & 7))) & 1) break;
            }
            vx = px*matr[0][0] + py*matr[0][1] + pz*matr[0][2];
            vz = px*matr[2][0] + py*matr[2][1] + pz*matr[2][2];
            v = (vx*vx + vz*vz);
            if((v <= rr2) && (v >= rr1)) {
               sx = px*shade[0][0] + py*shade[0][1] + pz*shade[0][2];
               sy = px*shade[1][0] + py*shade[1][1] + pz*shade[1][2];
               sz = px*shade[2][0] + py*shade[2][1] + pz*shade[2][2];
               s = (sx*sx + sy*sy);

               if((sz < 0) && (s <= rr))
                  calp = 0.0;
               else {
                  calp = (float)shad/100;
/* added 122697 to transition into shaded part */
                  if(sz < 0) {
                     dr = (float)sqrt((double)s)-(float)r;
                     if (dr <= (float)pa) calp *= dr/(float)pa;
                  }
/* */
               }
/* added 010198 to transition into planet atmosphere */
               if(sz < 0) {
                  dr = (float)r-(float)sqrt((double)(px*px + py*py));
                  if ((dr > 0) && (dr <= (float)pa)) calp *= ((float)pa-dr)/pa;
               }
/* */

               f = (float)sqrt((double)v);
               alp = (float)atan2((double)vz,(double)vx);

/* get color from map and place on ring plane */

               ringmap((f-(float)r1)/(r2-r1),alp,calp);
            }
            break;

            case PLNT:
            case MOON:

/* generate planet or moon (only the parts we can see though) */

            if(mcz < 0) {
	      if((abuf[bpos>>3] >> (7-((x+xoff) & 7))) & 1) break;
            }

            px = (x - (long)pcx)/aspect;

/* calculate for coordinates on planet surface only */

            pr2 = px*px + py*py;
            if((pz2 = rr2 - pr2) >= 0L) {
               pz = (float)sqrt((double)pz2);
               pr = (float)sqrt((double)pr2);

/* compute cosine of angle between normal and light source */

               vx = px*matr[0][0] + py*matr[0][1] + pz*matr[0][2];
               vy = px*matr[1][0] + py*matr[1][1] + pz*matr[1][2];
               vz = px*matr[2][0] + py*matr[2][1] + pz*matr[2][2];

               dr = (float)r2-pr;
               if (dr <= (float)pa) va = -255.0*((float)pa-dr)/pa;
               else va = tran;

               calp=(px*light[0][2]+py*light[1][2]+pz*light[2][2])/(float)r2;

               if(calp > 1.0) calp = 1.0;
               if(calp < 0.0) calp = 0.0;

/* compute perspective coefficient */

/*
               f = (float)(zres-pz)/(float)(pdist-mcz);
*/
               f = 1.0;

/* get color from planet map */

               map(vx,vy,vz,va,calp,f);
               if(!pixel) {
		 if((abuf[bpos>>3] >> (7-((x+xoff) & 7))) & 1) break;

/* compute cosine of angle between normal and light source */

                  vx = px*matr[0][0] + py*matr[0][1] - pz*matr[0][2];
                  vy = px*matr[1][0] + py*matr[1][1] - pz*matr[1][2];
                  vz = px*matr[2][0] + py*matr[2][1] - pz*matr[2][2];

                  calp=-(px*light[0][2]+py*light[1][2]+pz*light[2][2])/(float)r2;

                  if(calp > 1.0) calp = 1.0;
                  if(calp < 0.0) calp = 0.0;

/* compute perspective coefficient */

/*
                  f = (float)(zres+pz)/(float)(pdist-mcz);
*/
                  f = 1.0;

                  map(vx,vy,vz,va,calp,f);
               }
            }
            break;
         }
         if(pixel) {
	   if(!(type & PICT) && (dr > pa)) abuf[bpos>>3] |= (UBYTE)(1 << (7-((x+xoff) & 7)));
            rbuf[bpos] = rc;
            gbuf[bpos] = gc;
            bbuf[bpos] = bc;
         }
      }
   }

   made |= type;

   return(1);

} /* end generate */

void reset()
{
   identity(moon);
   identity(trans);
   identity(light);
   identity(shade);
   xres=zres=wres=vres=320L;
   yres=200L;
   xcen=0L;
   ycen=0L;
   aspect = 1.0;
   dang=0;
   pdist = 320;
   patmo = 0;
   mdist = 0;
   rad = prad = mrad = rad1 = rad2 = 100;
   shad = 50;
   ra = ga = ba = 0.0;
   rb = gb = bb = 0;
   wrap = 1;
   olap = 0;
   tran = 0;
   type = PLNT;
}

void close_map()
{
   if(mask & RMAP) {
      free(rmap);
      mask ^= RMAP;
   }
   if(mask & GMAP) {
      free(gmap);
      mask ^= GMAP;
   }
   if(mask & BMAP) {
      free(bmap);
      mask ^= BMAP;
   }
}

char pgetc(FILE *file)
{
    int ich;
    char ch;

    ich = (int)fgetc( file );
    if ( ich == EOF ) return(0);
    ch = (char) ich;
    
    if ( ch == '#' ) {
        do {
            ich = (int)fgetc( file );
            if ( ich == EOF ) return(0);
            ch = (char) ich;
        }
        while ( ch != '\n' );
    }

    return ch;
}

int getint(FILE *file)
{
    char ch;
    int i;

    do {
      ch = pgetc( file );
    } while ( ch == ' ' || ch == '\t' || ch == '\n' );

    if ( ch < '0' || ch > '9' ) return(0);

    i = 0;
    do {
        i = i * 10 + ch - '0';
        ch = pgetc( file );
    } while ( ch >= '0' && ch <= '9' );

    return i;
}

/* Open stream. If file specification begins with '|' then open as a piped
   stream otherwise open as a file */

FILE *Open(char *fspec, char *mode)
{
    if (*fspec == '|') {
        ++fspec;
	pipe = TRUE;
	return popen(fspec, mode);
    }
    else {
        pipe = FALSE;
	return fopen(fspec, mode);
    }
}

/* Close stream. */

int Close(FILE *stream)
{
    int rc = 0;

    if (pipe) rc = pclose(stream);
    else rc = fclose(stream);
    pipe = FALSE;

    return rc;
}

int ReadPPM(char *fspec)
{
    FILE *fp;
    int i,j,n,pos,Maxval;
    long ich1,ich2,id,buflen;

    if ((fp = Open(fspec,"r")) == NULL) {
        fprintf(stderr,"read - could not open map!!!\n");
        return(0);
    }

    ich1 = (long)fgetc(fp);
    if ( ich1 == EOF ) {
        fprintf(stderr,"read - premature EOF reading magic number\n");
	Close(fp);
        return(0);
    }
    ich2 = (long)fgetc(fp);
    if ( ich2 == EOF ) {
        fprintf(stderr,"read - premature EOF reading magic number\n");
	Close(fp);
        return(0);
    }
    id = ich1 * 256L + ich2;
    if ( (id != PGM_FORMAT) && (id != RPGM_FORMAT)
      && (id != PPM_FORMAT) && (id != RPPM_FORMAT)) {
        fprintf(stderr,"read - not a PPM or PGM file!!!\n");
	Close(fp);
        return(0);
    }

    Width = (int)getint(fp);
    Height = (int)getint(fp);

    Maxval = getint(fp);
    i = 256/(Maxval+1);
    n = 0;
    while(i > 1) { n++; i >>= 1; }

    buflen = (long)Width*Height;
    if ((rmap = (UBYTE *)malloc((size_t)buflen)) == NULL) {
        fprintf(stderr,"read - insufficient memory for map!!!\n");
	Close(fp);
        return(0);
    }
    mask |= RMAP;
    if ((gmap = (UBYTE *)malloc((size_t)buflen)) == NULL) {
        fprintf(stderr,"read - insufficient memory for map!!!\n");
        Close(fp);
        return(0);
    }
    mask |= GMAP;
    if ((bmap = (UBYTE *)malloc((size_t)buflen)) == NULL) {
        fprintf(stderr,"read - insufficient memory for map!!!\n");
	Close(fp);
        return(0);
    }
    mask |= BMAP;

    pos = 0;
    for (i=0;i<Height;i++) { /* process n lines/screen */
        switch(id) {
            case PGM_FORMAT:
                for (j=0;j<Width;j++) {
                    rmap[pos] = gmap[pos] = bmap[pos] = (UBYTE)(getint(fp) << n);
                    pos++;
                }
                break;
            case RPGM_FORMAT:
                for (j=0;j<Width;j++) {
                    rmap[pos] = gmap[pos] = bmap[pos] = (UBYTE)(fgetc(fp) << n);
                    pos++;
                }
                break;
            case PPM_FORMAT:
                for (j=0;j<Width;j++) {
                    rmap[pos] = (UBYTE)(getint(fp) << n);
                    gmap[pos] = (UBYTE)(getint(fp) << n);
                    bmap[pos] = (UBYTE)(getint(fp) << n);
                    pos++;
                }
                break;
            case RPPM_FORMAT:
                for (j=0;j<Width;j++) {
                    rmap[pos] = (UBYTE)(fgetc(fp) << n);
                    gmap[pos] = (UBYTE)(fgetc(fp) << n);
                    bmap[pos] = (UBYTE)(fgetc(fp) << n);
                    pos++;
                }
                break;
        }
    }

    while((long)fgetc(fp) != EOF);

    Close(fp);
    return(1);
} /* end ReadPPM */

void SavePPM(char *savefile)
{
    FILE *ifp;
    long i,j,pos;

    if(made) {
        if ((ifp = Open(savefile,"w")) == NULL) {
            fprintf(stderr,"save - could not open file!!!\n");
            return;
        }

        fprintf(ifp,"P6\n%ld %ld\n255\n",vres,yres);
        pos = 0;
        for(j=0;j<yres;j++) {
            for(i=0;i<vres;i++) {
                putc(rbuf[pos],ifp);
                putc(gbuf[pos],ifp);
                putc(bbuf[pos],ifp);
                pos++;
            }
        }
        fflush(ifp);
        Close(ifp);
    }
    else fprintf(stderr,"save - Nothing to save!!!\n");
}

#define MAXARG 16

APIRET APIENTRY planetHandler( 
   PRXSTRING command, 
   PUSHORT flags,
   PRXSTRING returnstring)
{
  char *args[MAXARG],*args0;
  char field[80];
  char *icom,*ocom;
  char commands[13][10] = { "generate","save","planet","moon","ring",
                            "light","map","image","center",
                            "reset","view","stop","status" };
  int ncom=13;
  /*  long primary=0; */
  int i,m,n,nc,argn;
  int r,g,b,ang;

  if (command->strptr != NULL) {
    for(i=0;i<MAXARG;i++) args[i] = NULL;

    args0 = (char *)malloc((size_t)(command->strlength + 1));
    if (args0 == NULL) { 
        fprintf(stderr,"main - could not allocate memory for args\n");
        return 0;
    }
    strcpy(args0, command->strptr);
    argn = parse(args0,args,MAXARG);

    m = 0;
    nc = ncom;
    for(i=0;i<ncom;i++) {
        icom = args[0];
        ocom = commands[i];
        n = 0;
        while(((*icom | ' ') == *ocom) && (*icom != 0) && (*ocom != 0)) {
            n++; icom++; ocom++;
        }
        if(*icom == 0) {
            if(n > m) { m = n; nc = i; }
            else if(n == m) nc = -1;
        }
    }

    field[0] = 0;

    switch(nc) {

/*******************************************************************************
   generate(depth) - generate object
   where:   depth = out-of-plane depth of object (only applies for 3D, see view)
   result:  none
*******************************************************************************/
        case 0:
            i = 0;
	    if(!open_things()) {
	      close_things();
	      fprintf(stderr,"main - could not allocate memory for bitmap\n");
	      return(0);
	    }
            if(argn > 1) sscanf(args[1],"%d",&i);
	    generate(0,i);
	    if (dang) generate(1,i);
        break;

/*******************************************************************************
   save(filename) - save image
   where:   filename = the name of the file to store image in
   result:  none
*******************************************************************************/
        case 1: SavePPM(args[1]); break;

/*******************************************************************************
   planet(radius,distance,atmo,angle(s),..) - create planet
   where:   radius   = the radius of the planet
            distance = the distance from viewpoint
            atmo     = thickness of atmosphere
            angle(s) = angles (degrees) preceded by axis of rotation
                       (x,y or z)
   result:  none
*******************************************************************************/
        case 2:
            pdist = wres;
            patmo = 0;
            mdist = 0;
            prad = mrad = 100;
            if(argn > 1) {
                if (sscanf(args[1],"%d",&i)) prad = i;
            }
            if(argn > 2) {
                if (sscanf(args[2],"%d",&i)) pdist = i;
            }
            if(argn > 3) {
                if (sscanf(args[3],"%d",&i)) patmo = i;
            }
            identity(trans);
            identity(moon);
            if(argn > 4) {
                for(i = 4; i < argn; ++i) {
                    if(sscanf(args[i],"%*c %d",&ang) == 1)
                        rotate(args[i][0] | ' ',-ang,trans);
                }
                for(i = argn-1; i > 3; i--) {
                    if(sscanf(args[i],"%*c %d",&ang) == 1)
                        rotate(args[i][0] | ' ',ang,moon);
                }
            }
            type = PLNT;
            rad = prad;
        break;

/*******************************************************************************
   moon(radius,distance,angle(s),..) - create moon
   where:   radius   = the radius of the moon
            distance = the distance from the planet
            angle(s) = upto three angles (degrees) preceded by axis of rotation
                       (x,y or z)
   result:  none
*******************************************************************************/
        case 3:
            patmo = 0;
            mdist = 0;
            mrad = 100;
            if(argn > 1) {
                if (sscanf(args[1],"%d",&i)) mrad = i;
            }
            if(argn > 2) {
                if (sscanf(args[2],"%d",&i)) mdist = i;
            }
            identity(trans);
            if(argn > 3) {
                for(i = 3; i < argn; ++i) {
                    if(sscanf(args[i],"%*c %d",&ang) == 1)
                        rotate(args[i][0] | ' ',-ang,trans);
                }
            }
            type = MOON;
            rad = mrad;
        break;

/*******************************************************************************
   ring(radius1,radius2,shade) - create ring
   where:   radius1 = the inner radius of the ring
            radius2 = the outer radius of the ring
            shade   = the shade factor (0-100)
   result:  none
*******************************************************************************/
        case 4:
            rad1 = rad2 = rad;
            shad = 50;
            if(argn > 1) {
                if (sscanf(args[1],"%d",&i)) rad1 = i;
            }
            if(argn > 2) {
                if (sscanf(args[2],"%d",&i)) rad2 = i;
            }
            if(argn > 3) sscanf(args[3],"%d",&shad);
            type = RING;
            if(rad1 < rad) rad1 = rad;
            if(rad2 < rad1) rad2 = rad1;
        break;

/*******************************************************************************
   light(diffuse,angle(s),...) - position light source
   where:   diffuse = six digit hex value for color of diffuse light (rrggbb)
            angle(s) = upto three angles (degrees) preceded by axis of rotation
                       (x,y or z)
   result:  none
*******************************************************************************/
        case 5:
            ra = ga = ba = 0.0;
            ang = 0;
            if(argn > 1) {
                i = sscanf(args[1],"%2x%2x%2x",&r,&g,&b);
                if(i > 0) ra = (float)r/255.0;
                else { r = 0; ra = 0.0; }
                if(i > 1) ga = (float)g/255.0;
                else { g = 0; ga = 0.0; }
                if(i > 2) ba = (float)b/255.0;
                else { b = 0; ba = 0.0; }
            }
            identity(light);
            identity(shade);
            if(argn > 2) {
                for(i = 2; i < argn; ++i) {
                    if(sscanf(args[i],"%*c %d",&ang) == 1) {
                        rotate(args[i][0] | ' ',-ang,shade);
                    }
                }
                for(i = argn-1; i > 1; i--) {
                    if(sscanf(args[i],"%*c %d",&ang) == 1) {
                        rotate(args[i][0] | ' ',ang,light);
                    }
                }
            }
        break;

/*******************************************************************************
   map(filename,wrap,olap,transparency,threshold) - specify map
   where:   filename     = the name of the file for the map
            wrap         = 'full' or 'half' wrap
            olap         = the amount of overlap at seam
            transparency = the degree of transparency of the map (0 - 255)
            threshold    = six digit hex value for color threshold (rrggbb)
   result:  none
*******************************************************************************/
        case 6:
            rb = gb = bb = 0;
            wrap = 2;
            olap = 0;
            tran = 0;
            if((argn > 2) && ((args[2][0] | ' ') == 'f')) wrap = 2;
            else wrap = 1;
            if(argn > 3) sscanf(args[3],"%d",&olap);
            if(argn > 4) sscanf(args[4],"%d",&tran);
            if(argn > 5) sscanf(args[5],"%2x%2x%2x",&rb,&gb,&bb);
            close_map();
            if(!ReadPPM(args[1])) {
                close_map();
                fprintf(stderr,"main - error reading map!!!\n");
            }
        break;

/*******************************************************************************
   image - specify image
   result:  none
*******************************************************************************/
        case 7: type = PICT; break;

/*******************************************************************************
   center(x,y) - center planet or image about x,y (0,0 is center of screen)
   where:   x = x coordinate of center
            y = y coordinate of center
   result:  none
*******************************************************************************/
        case 8:
            if(argn > 1) {
                if (sscanf(args[1],"%d",&i)) xcen = i;
            }
            else xcen = 0;
            if(argn > 2) {
                if (sscanf(args[2],"%d",&i)) ycen = i;
            }
            else ycen = 0;
        break;

/*******************************************************************************
   reset(wres) - reset all parameters and free up memory
   where:   wres    = width in world coord
   result:  none
*******************************************************************************/
        case 9:
            close_map();
            close_things();
            made = 0;
            reset();
            if(argn > 1) sscanf(args[1],"%ld",&wres);
        break;

/*******************************************************************************
   view(xres,yres,dang,aspect) - set view
   where:   xres    = x resolution of view
            yres    = y resolution of view
            dang    = angle shift in view (positive non-zero value implies 3D)
            aspect  = aspect ratio of view
   result:  none
*******************************************************************************/
        case 10:
            close_things();
            made = 0;
            if(argn > 1) sscanf(args[1],"%ld",&xres);
            if(argn > 2) sscanf(args[2],"%ld",&yres);
            if(argn > 3) sscanf(args[3],"%d",&dang);
            if(argn > 4) sscanf(args[4],"%g",&aspect);
            zres = xres;
        break;

/*******************************************************************************
   stop - remove window and shut down
   result:  none
*******************************************************************************/
        case 11: close_down = TRUE; break;

        case 12: sprintf(field,"awaiting command"); break;

        case 13: fprintf(stderr,"main - unknown command\n"); break;

        case -1: fprintf(stderr,"main - ambiguous command\n"); break;
    }

    if (field[0]) {
        MAKERXSTRING(result.shvvalue,field,strlen(field));
        RexxVariablePool(&result);
    }
    returnstring->strptr = NULL;
    returnstring->strlength = 0;
    if (args0) free(args0);
  }
  return 0;
}

int main( int argc, char *argv[] )
{
  int rc ;
  short returnCode;
  RXSTRING Result ;

  if(argc > 1) {
    made = 0;
    reset();
    rc = RexxRegisterSubcomExe( hostname, (PFN) planetHandler, NULL ) ;

    Result.strlength = 200 ;
    Result.strptr = malloc( 200 ) ;

    rc = RexxStart( 0, NULL, argv[1], 0, hostname, RXCOMMAND,
                    NULL, &returnCode, &Result ) ;
    if (rc < 0) rc = -rc;
  }
  return 0;
}
