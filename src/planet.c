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

 Copyright � 1990,1991,1992 by Russell Leighton

 ==============================================================================*/

#include "planet.h"

extern int parse(char *input, char *args[], int narg);

double vabs(double value) {
    return value < 0.0 ? -value : value;
}

/* Equate two matrices M = N */

void mset(float (*n)[3], float (*m)[3]) {
    int i, j;

    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++)
            m[i][j] = n[i][j];
    }
}

/* Multiply two matrices together M = N*M */

void mmul(float (*n)[3], float (*m)[3]) {
    int i, j, k;
    float p[3][3];

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            p[i][j] = 0.0;
            for (k = 0; k < 3; k++)
                p[i][j] += n[i][k] * m[k][j];
        }
    }
    mset(p, m);
}

/* Divide two matrices together M = M/N */

void mdiv(float (*n)[3], float (*m)[3]) {
    int i, j;
    float p[3][3], det = 0.0;

    for (i = 0; i < 3; i++) {
        det += (n[0][i]
                * (n[1][(i + 1) % 3] * n[2][(i + 2) % 3]
                        - n[1][(i + 2) % 3] * n[2][(i + 1) % 3]));
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            p[i][j] =
                    ((n[(i + 1) % 3][(j + 1) % 3] * n[(i + 2) % 3][(j + 2) % 3])
                            - (n[(i + 1) % 3][(j + 2) % 3]
                                    * n[(i + 2) % 3][(j + 1) % 3])) / det;
        }
    }
    mmul(p, m);
}

/* Equate two vectors M = N */

void vset(float n[3], float m[3]) {
    int i;

    for (i = 0; i < 3; i++)
        m[i] = n[i];
}

/* Multiply a vector and matrix together M = N*M */

void vmul(float (*n)[3], float m[3]) {
    int i, j;
    float p[3];

    for (i = 0; i < 3; i++) {
        p[i] = 0.0;
        for (j = 0; j < 3; j++) {
            p[i] += n[i][j] * m[j];
        }
    }
    vset(p, m);
}

/* compute transformation and concatenate with previous transformation(s) */

void rotate(int axis, int ang, float (*m)[3]) {
    int i, j, k;
    float n[3][3], a;

    if (!ang)
        return;

    /* setup transformation matrix based on selected axis */

    switch (axis) {
    case 'x':
        i = 0;
        j = 1;
        k = 2;
        break;
    case 'y':
        i = 1;
        j = 2;
        k = 0;
        break;
    case 'z':
        i = 2;
        j = 0;
        k = 1;
        break;
    default:
        i = 0;
        j = 1;
        k = 2;
        break;
    }

    a = ang > 0 ? (float) ang + 0.5 : (float) ang - 0.5;
    n[i][i] = 1.0;
    n[j][j] = n[k][k] = (float) cos((double) a * PI / 180.0);
    n[i][j] = n[i][k] = n[j][i] = n[k][i] = 0.0;
    n[k][j] = (float) sin((double) a * PI / 180.0);
    n[j][k] = -n[k][j];

    mmul(n, m);

} /* end rotate */

/* Return identity matrix M = I */

void identity(float (*m)[3]) {
    int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i == j)
                m[i][j] = 1;
            else
                m[i][j] = 0;
        }
    }
}

/* get color from map at specified px,py,pz */

void map(float v[3], int da, float ca) {
    long x, y;
    long pos;
    int r1, g1, b1, a1, r2, g2, b2;
    float lon, lat, cf;

    /* compute longitude and latitude */
    lon = (float) atan2((double) v[0], (double) v[2]);
    if (wrap == 1) {
        lon *= 2.0;
        if (lon < 0.0)
            lon += PI;
        else
            lon -= PI;
    }
    lat = (float) atan2((double) v[1], sqrt((double) (v[0] * v[0] + v[2] * v[2])));

    x = (.5 + lon / (2.0 * PI)) * (Width - olap - 1);
    y = (.5 - lat / PI) * (Height - 1);

    pos = y * Width + x;
    r1 = (int) rmap[pos];
    g1 = (int) gmap[pos];
    b1 = (int) bmap[pos];
    a1 = da >= 0 ? da : (int) amap[pos];
    if ((x < olap) && (olap != 0)) {
        pos += Width - olap;
        r2 = (int) rmap[pos];
        g2 = (int) gmap[pos];
        b2 = (int) bmap[pos];
        r1 = (x * r1 + (olap - x) * r2) / olap;
        g1 = (x * g1 + (olap - x) * g2) / olap;
        b1 = (x * b1 + (olap - x) * b2) / olap;
    }

    if (ca > ra)
        r2 = ca * r1;
    else
        r2 = ra * r1;
    if (ca > ga)
        g2 = ca * g1;
    else
        g2 = ga * g1;
    if (ca > ba)
        b2 = ca * b1;
    else
        b2 = ba * b1;
    cf = a1 / 255.0;
    rc = cf * r2 + (1.0 - cf) * rc;
    gc = cf * g2 + (1.0 - cf) * gc;
    bc = cf * b2 + (1.0 - cf) * bc;
    pixel = TRUE;
} /* end map */

void ringmap(float r, float ang, float ca, float af, int shadow) {
    long x, y;
    long pos;
    int r1, g1, b1, r2 = 0, g2 = 0, b2 = 0;
    float cf;

    x = (ang / PI + 1.0) * (Width - 1) / 2.0;
    y = r * (Height - 1);

    pos = y * Width + x;
    if (!shadow) {
        r1 = (int) rmap[pos];
        g1 = (int) gmap[pos];
        b1 = (int) bmap[pos];
        if (x < olap) {
            pos += Width - olap;
            r2 = (int) rmap[pos];
            g2 = (int) gmap[pos];
            b2 = (int) bmap[pos];
            r1 = (x * r1 + (olap - x) * r2) / olap;
            g1 = (x * g1 + (olap - x) * g2) / olap;
            b1 = (x * b1 + (olap - x) * b2) / olap;
        }

        if (ca > ra)
            r2 = ca * r1;
        else
            r2 = ra * r1;
        if (ca > ga)
            g2 = ca * g1;
        else
            g2 = ga * g1;
        if (ca > ba)
            b2 = ca * b1;
        else
            b2 = ba * b1;
    }
    cf = af * (float) amap[pos] / 255.0;
    rc = cf * r2 + (1.0 - cf) * rc;
    gc = cf * g2 + (1.0 - cf) * gc;
    bc = cf * b2 + (1.0 - cf) * bc;
    pixel = TRUE;
} /* end ringmap */

int open_things() {
    long buflen;

    vres = xres;
    if (dang)
        vres = 2 * xres;
    buflen = vres * yres;

    if (!(mask & RED)) {
        if ((rbuf = calloc((long) buflen, 1L)) == NULL) {
            fprintf(stderr, "main - insufficient memory!!!\n");
            return (0);
        }
        mask |= RED;
    }
    if (!(mask & GREEN)) {
        if ((gbuf = calloc((long) buflen, 1L)) == NULL) {
            fprintf(stderr, "main - insufficient memory!!!\n");
            return (0);
        }
        mask |= GREEN;
    }
    if (!(mask & BLUE)) {
        if ((bbuf = calloc((long) buflen, 1L)) == NULL) {
            fprintf(stderr, "main - insufficient memory!!!\n");
            return (0);
        }
        mask |= BLUE;
    }
    buflen = vres / 8L * yres;
    if (!(mask & ALPHA)) {
        if ((abuf = calloc((long) buflen, 1L)) == NULL) {
            fprintf(stderr, "main - insufficient memory!!!\n");
            return (0);
        }
        mask |= ALPHA;
    }
    return (1);
}

void close_things() {
    if (mask & RED) {
        free(rbuf);
        mask ^= RED;
    }
    if (mask & GREEN) {
        free(gbuf);
        mask ^= GREEN;
    }
    if (mask & BLUE) {
        free(bbuf);
        mask ^= BLUE;
    }
    if (mask & ALPHA) {
        free(abuf);
        mask ^= ALPHA;
    }
}

int generate(int lflag, int depth) {
    int d, da, r = 0, r1 = 0, r2 = 0;
    long pcx = 0, pcy = 0;
    long x, y, bpos, mpos;
    long s2, v2, pz2, rr = 0, rr1 = 0, rr2 = 0, pr2;
    long mcx, mcy, mcz = 0;
    long minx, maxx, miny, maxy, xoff;
    float f, calp, fdist, alp, pa, pr, dr, af;
    float v[3], p[3], s[3];
    float matr[3][3], ring[3][3];

    /* determine minimum possible window to conserve computation */
    if (type & PICT) {
        minx = 0L;
        miny = 0L;
        maxx = xres;
        maxy = yres;
    } else {
        if (type & MOON) {
            mset(moon, matr);
            if (dang) {
                if (lflag) rotate('y', -(dang >> 1), matr);
                else rotate('y', (dang >> 1), matr);
            }
            s2 = (mdist * xres) / pdist;
            mcx = s2 * matr[0][2] * aspect;
            mcy = s2 * matr[1][2];
            mcz = s2 * matr[2][2];
        } else
            mcx = mcy = mcz = 0;
        d = (pdist * xres) / pdist - mcz;
        if (d <= 0)
            return (0);
        fdist = (float) zres / (float) d;
        r = (rad * xres * fdist) / pdist;
        pcx = (xcen * xres) / pdist + xres / 2;
        if (dang) {
            if (lflag) pcx -= depth >> 1;
            else pcx += depth >> 1;
        }
        pcy = (ycen * xres) / pdist + yres / 2;
        pcx += mcx;
        pcy -= mcy;
        if (type & RING) {
            r1 = (rad1 * xres * fdist) / pdist;
            r2 = (rad2 * xres * fdist) / pdist;
            rr1 = (long) r1 * r1;
            rr = r * r;
        } else
            r2 = r;
        rr2 = (long) r2 * r2;

        minx = (long) pcx - (long) (r2 * aspect);
        minx = (minx > 0L) ? minx : 0L;
        miny = (long) pcy - (long) r2;
        miny = (miny > 0L) ? miny : 0L;
        maxx = (long) pcx + (long) (r2 * aspect);
        maxx = (maxx < xres) ? maxx : xres;
        maxy = (long) pcy + (long) r2;
        maxy = (maxy < yres) ? maxy : yres;
        identity(matr);
        if (dang) {
            if (lflag) rotate('y', (dang >> 1), matr);
            else rotate('y', -(dang >> 1), matr);
        }
        mmul(trans, matr);
    }
    pa = (patmo * xres) / pdist;

    if (lflag)
        xoff = xres;
    else
        xoff = 0;

    /* loop on destination display coordinates */
    for (y = miny; y < maxy; y++) {
        if (type & PICT) {
            p[1] = (y * (long) Height) / yres;
        } else
            p[1] = (long) pcy - y;

        for (x = minx; x < maxx; ++x) {

            pixel = FALSE;
            bpos = vres * y + (x + xoff);
            rc = rbuf[bpos];
            gc = gbuf[bpos];
            bc = bbuf[bpos];
            dr = pa;

            switch (type) {
            case PICT:
                p[0] = x;
                if (dang) {
                    if (lflag) p[0] += depth >> 1;
                    else p[0] -= depth >> 1;
                }
                p[0] *= (float) Width / xres;
                if ((p[0] < 0) || (p[0] >= Width)) {
                    rc = gc = bc = 0;
                    break;
                }
                mpos = p[1] * Width + p[0];
                rc = rmap[mpos];
                gc = gmap[mpos];
                bc = bmap[mpos];
                pixel = TRUE;
                break;

            case RING:
                p[0] = (x - (long) pcx) / aspect;
                pr2 = p[0] * p[0] + p[1] * p[1];

                mset(matr,ring);
                mdiv(light,ring);

                // render ring shadow on planet surface
                if ((pz2 = rr - pr2) >= 0L && vabs(ring[2][1]) > 0.0001) {
                    p[2] = (float) sqrt((double) pz2);
                    vset(p,v);
                    vmul(matr,v);
                    v2 = (v[0] * v[0] + v[2] * v[2]);
                    f = (float) sqrt((double) v2) + (v[1] * ring[2][2] / ring[2][1]);

                    if ((f >= r1) && (f <= r2)) {
                        calp = (p[0] * light[0][2] + p[1] * light[1][2]
                                + p[2] * light[2][2]) / (float) r;
                        if (calp > 1.0) calp = 1.0;
                        if (calp < 0.0) calp = 0.0;
                        alp = (float) atan2((double) v[2], (double) v[0]);
                        af = calp * ((float) shad / 100);
                        ringmap((f - (float) r1) / (r2 - r1), alp, calp, af, TRUE);
                    }
                }
                if (vabs(matr[1][2]) <= 0.0001)
                    break;
                p[2] = -(p[0] * matr[1][0] + p[1] * matr[1][1]) / matr[1][2];
                vset(p,v);
                vmul(matr,v);
                v2 = (v[0] * v[0] + v[2] * v[2]);
                if ((v2 <= rr2) && (v2 >= rr1) &&
                    !((p[2] < 0) && ((abuf[bpos >> 3] >> (7 - ((x + xoff) & 7))) & 1))) {
                    vset(p,s);
                    vmul(shade,s);
                    s2 = (s[0] * s[0] + s[1] * s[1]);

                    if ((s[2] < 0) && (s2 <= rr))
                        calp = 0.0;
                    else {
                        calp = (float) shad / 100;
                        /* transition into shaded part */
                        if (s[2] < 0) {
                            dr = (float) sqrt((double) s2) - (float) r;
                            if (dr <= pa)
                                calp *= dr / pa;
                        }
                    }
                    /* transition into planet atmosphere */
                    af = 1.0;
                    if (p[2] < 0) {
                        dr = (float) r - (float) sqrt((double) (p[0] * p[0] + p[1] * p[1]));
                        if ((dr > 0) && (dr <= pa))
                            af = (pa - dr) / pa;
                    }

                    f = (float) sqrt((double) v2);
                    alp = (float) atan2((double) v[2], (double) v[0]);

                    /* get color from map and place on ring plane */
                    ringmap((f - (float) r1) / (r2 - r1), alp, calp, af, FALSE);
                }
                break;

            case PLNT:
            case MOON:

                /* generate planet or moon (only the parts we can see though) */
                if (mcz < 0) {
                    if ((abuf[bpos >> 3] >> (7 - ((x + xoff) & 7))) & 1)
                        break;
                }

                p[0] = (x - (long) pcx) / aspect;

                /* calculate for coordinates on planet surface only */
                pr2 = p[0] * p[0] + p[1] * p[1];
                if ((pz2 = rr2 - pr2) >= 0L) {
                    p[2] = (float) sqrt((double) pz2);
                    pr = (float) sqrt((double) pr2);

                    vset(p,v);
                    vmul(matr,v);
                    dr = (float) r2 - pr;

                    if (dr <= pa)
                        da = 255.0 * (1.0 - (pa - dr) / pa);
                    else
                        da = -1.0;

                    /* compute cosine of angle between normal and light source */
                    calp = (p[0] * light[0][2] + p[1] * light[1][2]
                            + p[2] * light[2][2]) / (float) r2;

                    if (calp > 1.0) calp = 1.0;
                    if (calp < 0.0) calp = 0.0;

                    /* get color from planet map */
                    map(v, da, calp);
                }
                break;
            }
            if (pixel) {
                if (!(type & PICT) && (dr > pa))
                    abuf[bpos >> 3] |= (UBYTE) (1 << (7 - ((x + xoff) & 7)));
                rbuf[bpos] = rc;
                gbuf[bpos] = gc;
                bbuf[bpos] = bc;
            }
        }
    }

    made |= type;

    return (1);

} /* end generate */

void reset() {
    identity(moon);
    identity(trans);
    identity(light);
    identity(shade);
    xres = zres = vres = 320L;
    yres = 200L;
    xcen = 0L;
    ycen = 0L;
    aspect = 1.0;
    dang = 0;
    pdist = 320;
    patmo = 0;
    mdist = 0;
    rad = prad = mrad = rad1 = rad2 = 100;
    shad = 50;
    ra = ga = ba = 0.0;
    wrap = 1;
    olap = 0;
    type = PLNT;
}

void close_map() {
    if (mask & RMAP) {
        free(rmap);
        mask ^= RMAP;
    }
    if (mask & GMAP) {
        free(gmap);
        mask ^= GMAP;
    }
    if (mask & BMAP) {
        free(bmap);
        mask ^= BMAP;
    }
    if (mask & AMAP) {
        free(amap);
        mask ^= AMAP;
    }
}

char pgetc(FILE *file) {
    int ich;
    char ch;

    ich = (int) fgetc(file);
    if (ich == EOF)
        return (0);
    ch = (char) ich;

    if (ch == '#') {
        do {
            ich = (int) fgetc(file);
            if (ich == EOF)
                return (0);
            ch = (char) ich;
        } while (ch != '\n');
    }

    return ch;
}

int getint(FILE *file) {
    char ch;
    int i;

    do {
        ch = pgetc(file);
    } while (ch == ' ' || ch == '\t' || ch == '\n');

    if (ch < '0' || ch > '9')
        return (0);

    i = 0;
    do {
        i = i * 10 + ch - '0';
        ch = pgetc(file);
    } while (ch >= '0' && ch <= '9');

    return i;
}

/* Open stream. If file specification begins with '|' then open as a piped
 stream otherwise open as a file */

FILE *Open(char *fspec, char *mode) {
    if (*fspec == '|') {
        ++fspec;
        pipe = TRUE;
        return popen(fspec, mode);
    } else {
        pipe = FALSE;
        return fopen(fspec, mode);
    }
}

/* Close stream. */

int Close(FILE *stream) {
    int rc = 0;

    if (pipe)
        rc = pclose(stream);
    else
        rc = fclose(stream);
    pipe = FALSE;

    return rc;
}

int ReadMap(char *fspec) {
    FILE *fp;
    int i, j, n, pos, Maxval;
    long ich1, ich2, id, buflen;

    if ((fp = Open(fspec, "r")) == NULL) {
        fprintf(stderr, "read - could not open map!!!\n");
        return (0);
    }

    ich1 = (long) fgetc(fp);
    if (ich1 == EOF) {
        fprintf(stderr, "read - premature EOF reading magic number\n");
        Close(fp);
        return (0);
    }
    ich2 = (long) fgetc(fp);
    if (ich2 == EOF) {
        fprintf(stderr, "read - premature EOF reading magic number\n");
        Close(fp);
        return (0);
    }
    id = ich1 * 256L + ich2;
    if ((id != PGM_FORMAT) && (id != RPGM_FORMAT) && (id != PPM_FORMAT)
            && (id != RPPM_FORMAT)) {
        fprintf(stderr, "read - not a PPM or PGM file!!!\n");
        Close(fp);
        return (0);
    }

    Width = (int) getint(fp);
    Height = (int) getint(fp);

    Maxval = getint(fp);
    i = 256 / (Maxval + 1);
    n = 0;
    while (i > 1) {
        n++;
        i >>= 1;
    }

    buflen = (long) Width * Height;
    if ((rmap = (UBYTE *) malloc((long) buflen)) == NULL) {
        fprintf(stderr, "read - insufficient memory for map!!!\n");
        Close(fp);
        return (0);
    }
    mask |= RMAP;
    if ((gmap = (UBYTE *) malloc((long) buflen)) == NULL) {
        fprintf(stderr, "read - insufficient memory for map!!!\n");
        Close(fp);
        return (0);
    }
    mask |= GMAP;
    if ((bmap = (UBYTE *) malloc((long) buflen)) == NULL) {
        fprintf(stderr, "read - insufficient memory for map!!!\n");
        Close(fp);
        return (0);
    }
    mask |= BMAP;
    if ((amap = (UBYTE *) malloc((long) buflen)) == NULL) {
        fprintf(stderr, "read - insufficient memory for map!!!\n");
        Close(fp);
        return (0);
    }
    mask |= AMAP;

    pos = 0;
    for (i = 0; i < Height; i++) { /* process n lines/screen */
        switch (id) {
        case PGM_FORMAT:
            for (j = 0; j < Width; j++) {
                rmap[pos] = gmap[pos] = bmap[pos] = (UBYTE) (getint(fp) << n);
                amap[pos] = 0xFF;
                pos++;
            }
            break;
        case RPGM_FORMAT:
            for (j = 0; j < Width; j++) {
                rmap[pos] = gmap[pos] = bmap[pos] = (UBYTE) (fgetc(fp) << n);
                amap[pos] = 0xFF;
                pos++;
            }
            break;
        case PPM_FORMAT:
            for (j = 0; j < Width; j++) {
                rmap[pos] = (UBYTE) (getint(fp) << n);
                gmap[pos] = (UBYTE) (getint(fp) << n);
                bmap[pos] = (UBYTE) (getint(fp) << n);
                amap[pos] = 0xFF;
                pos++;
            }
            break;
        case RPPM_FORMAT:
            for (j = 0; j < Width; j++) {
                rmap[pos] = (UBYTE) (fgetc(fp) << n);
                gmap[pos] = (UBYTE) (fgetc(fp) << n);
                bmap[pos] = (UBYTE) (fgetc(fp) << n);
                amap[pos] = 0xFF;
                pos++;
            }
            break;
        }
    }

    while ((long) fgetc(fp) != EOF)
        ;

    Close(fp);
    return (1);
} /* end ReadPPM */


int ReadAlphaMap(char *fspec) {
    FILE *fp;
    int i, j, n, w, h, pos, Maxval;
    long ich1, ich2, id;

    if (!(mask & AMAP)) {
        fprintf(stderr, "read - map must be defined before alpha map!!!\n");
        return (0);
    }

    if ((fp = Open(fspec, "r")) == NULL) {
        fprintf(stderr, "read - could not open map!!!\n");
        return (0);
    }

    ich1 = (long) fgetc(fp);
    if (ich1 == EOF) {
        fprintf(stderr, "read - premature EOF reading magic number\n");
        Close(fp);
        return (0);
    }
    ich2 = (long) fgetc(fp);
    if (ich2 == EOF) {
        fprintf(stderr, "read - premature EOF reading magic number\n");
        Close(fp);
        return (0);
    }
    id = ich1 * 256L + ich2;
    if ((id != PGM_FORMAT) && (id != RPGM_FORMAT)) {
        fprintf(stderr, "read - not a PGM file!!!\n");
        Close(fp);
        return (0);
    }

    w = (int) getint(fp);
    h = (int) getint(fp);

    if (w != Width || h != Height) {
        fprintf(stderr, "read - alpha map not same size as map!!!\n");
        Close(fp);
        return (0);
    }

    Maxval = getint(fp);
    i = 256 / (Maxval + 1);
    n = 0;
    while (i > 1) {
        n++;
        i >>= 1;
    }

    pos = 0;
    for (i = 0; i < Height; i++) { /* process n lines/screen */
        switch (id) {
        case PGM_FORMAT:
            for (j = 0; j < Width; j++) {
                amap[pos] = (UBYTE) (getint(fp) << n);
                pos++;
            }
            break;
        case RPGM_FORMAT:
            for (j = 0; j < Width; j++) {
                amap[pos] = (UBYTE) (fgetc(fp) << n);
                pos++;
            }
            break;
        }
    }

    while ((long) fgetc(fp) != EOF)
        ;

    Close(fp);
    return (1);
} /* end ReadAMAP */

void SavePPM(char *savefile) {
    FILE *ifp;
    long i, j, pos;

    if (made) {
        if ((ifp = Open(savefile, "w")) == NULL) {
            fprintf(stderr, "save - could not open file!!!\n");
            return;
        }

        fprintf(ifp, "P6\n%ld %ld\n255\n", vres, yres);
        pos = 0;
        for (j = 0; j < yres; j++) {
            for (i = 0; i < vres; i++) {
                putc(rbuf[pos], ifp);
                putc(gbuf[pos], ifp);
                putc(bbuf[pos], ifp);
                pos++;
            }
        }
        fflush(ifp);
        Close(ifp);
    } else
        fprintf(stderr, "save - Nothing to save!!!\n");
}

#define MAXARG 16

APIRET APIENTRY planetHandler(PRXSTRING command, PUSHORT flags,
        PRXSTRING returnstring) {
    char *args[MAXARG], *args0;
    char field[80];
    char *icom, *ocom;
    char commands[11][10] =
            { "generate", "save", "planet", "moon", "ring", "light", "map",
              "alpha", "image", "center", "view" };
    int ncom = 13;
    /*  long primary=0; */
    int i, m, n, nc, argn;
    int r, g, b, ang;

    if (command->strptr != NULL) {
        for (i = 0; i < MAXARG; i++)
            args[i] = NULL;

        args0 = (char *) malloc((long)(command->strlength + 1));
        if (args0 == NULL) {
            fprintf(stderr, "main - could not allocate memory for args\n");
            return 0;
        }
        strcpy(args0, command->strptr);
        argn = parse(args0, args, MAXARG);

        m = 0;
        nc = ncom;
        for (i = 0; i < ncom; i++) {
            icom = args[0];
            ocom = commands[i];
            n = 0;
            while (((*icom | ' ') == *ocom) && (*icom != 0) && (*ocom != 0)) {
                n++;
                icom++;
                ocom++;
            }
            if (*icom == 0) {
                if (n > m) {
                    m = n;
                    nc = i;
                } else if (n == m)
                    nc = -1;
            }
        }

        field[0] = 0;

        switch (nc) {

        /*******************************************************************************
         generate(depth) - generate object
         where:   depth = out-of-plane depth of object (only applies for 3D, see view)
         result:  none
         *******************************************************************************/
        case 0:
            i = 0;
            if (!open_things()) {
                close_things();
                fprintf(stderr,
                        "main - could not allocate memory for bitmap\n");
                return (0);
            }
            if (argn > 1)
                sscanf(args[1], "%d", &i);
            generate(0, i);
            if (dang)
                generate(1, i);
            break;

            /*******************************************************************************
             save(filename) - save image
             where:   filename = the name of the file to store image in
             result:  none
             *******************************************************************************/
        case 1:
            SavePPM(args[1]);
            break;

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
            pdist = 100;
            patmo = 0;
            mdist = 0;
            prad = mrad = 100;
            if (argn > 1) {
                if (sscanf(args[1], "%d", &i))
                    prad = i;
            }
            if (argn > 2) {
                if (sscanf(args[2], "%d", &i))
                    pdist = i;
            }
            if (argn > 3) {
                if (sscanf(args[3], "%d", &i))
                    patmo = i;
            }
            identity(trans);
            identity(moon);
            if (argn > 4) {
                for (i = 4; i < argn; ++i) {
                    if (sscanf(args[i], "%*c %d", &ang) == 1)
                        rotate(args[i][0] | ' ', -ang, trans);
                }
                for (i = argn - 1; i > 3; i--) {
                    if (sscanf(args[i], "%*c %d", &ang) == 1)
                        rotate(args[i][0] | ' ', ang, moon);
                }
            }
            type = PLNT;
            rad = prad;
            break;

            /*******************************************************************************
             moon(radius,distance,angle(s),..) - create moon
             where:   radius   = the radius of the moon
             distance = the distance from the planet
             atmo     = thickness of atmosphere
             angle(s) = angles (degrees) preceded by axis of rotation
             (x,y or z)
             result:  none
             *******************************************************************************/
        case 3:
            patmo = 0;
            mdist = 0;
            mrad = 100;
            if (argn > 1) {
                if (sscanf(args[1], "%d", &i))
                    mrad = i;
            }
            if (argn > 2) {
                if (sscanf(args[2], "%d", &i))
                    mdist = i;
            }
            if (argn > 3) {
                if (sscanf(args[3], "%d", &i))
                    patmo = i;
            }
            identity(trans);
            if (argn > 4) {
                for (i = 4; i < argn; ++i) {
                    if (sscanf(args[i], "%*c %d", &ang) == 1)
                        rotate(args[i][0] | ' ', -ang, trans);
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
            if (argn > 1) {
                if (sscanf(args[1], "%d", &i))
                    rad1 = i;
            }
            if (argn > 2) {
                if (sscanf(args[2], "%d", &i))
                    rad2 = i;
            }
            if (argn > 3)
                sscanf(args[3], "%d", &shad);
            type = RING;
            if (rad1 < rad)
                rad1 = rad;
            if (rad2 < rad1)
                rad2 = rad1;
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
            if (argn > 1) {
                i = sscanf(args[1], "%2x%2x%2x", &r, &g, &b);
                if (i > 0)
                    ra = (float) r / 255.0;
                else {
                    r = 0;
                    ra = 0.0;
                }
                if (i > 1)
                    ga = (float) g / 255.0;
                else {
                    g = 0;
                    ga = 0.0;
                }
                if (i > 2)
                    ba = (float) b / 255.0;
                else {
                    b = 0;
                    ba = 0.0;
                }
            }
            identity(light);
            identity(shade);
            if (argn > 2) {
                for (i = 2; i < argn; ++i) {
                    if (sscanf(args[i], "%*c %d", &ang) == 1) {
                        rotate(args[i][0] | ' ', -ang, shade);
                    }
                }
                for (i = argn - 1; i > 1; i--) {
                    if (sscanf(args[i], "%*c %d", &ang) == 1) {
                        rotate(args[i][0] | ' ', ang, light);
                    }
                }
            }
            break;

            /*******************************************************************************
             map(filename,wrap,olap) - specify map
             where:   filename     = the name of the file for the map
             wrap         = 'full' or 'half' wrap
             olap         = the amount of overlap at seam
             result:  none
             *******************************************************************************/
        case 6:
            wrap = 1;
            olap = 0;
            if ((argn > 2) && ((args[2][0] | ' ') == 'h'))
                wrap = 1;
            else
                wrap = 2;
            if (argn > 3)
                sscanf(args[3], "%d", &olap);
            close_map();
            if (!ReadMap(args[1])) {
                close_map();
                fprintf(stderr, "main - error reading map!!!\n");
            }
            break;

            /*******************************************************************************
             alpha(filename) - specify alpha map
             where:   filename     = the name of the file for the map
             result:  none
             *******************************************************************************/
        case 7:
            if (!ReadAlphaMap(args[1])) {
                close_map();
                fprintf(stderr, "main - error reading map!!!\n");
            }
            break;

            /*******************************************************************************
             image - specify image
             result:  none
             *******************************************************************************/
        case 8:
            type = PICT;
            break;

            /*******************************************************************************
             center(x,y) - center planet or image about x,y (0,0 is center of screen)
             where:   x = x coordinate of center
             y = y coordinate of center
             result:  none
             *******************************************************************************/
        case 9:
            if (argn > 1) {
                if (sscanf(args[1], "%d", &i))
                    xcen = i;
            } else
                xcen = 0;
            if (argn > 2) {
                if (sscanf(args[2], "%d", &i))
                    ycen = i;
            } else
                ycen = 0;
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
            close_map();
            close_things();
            made = 0;
            reset();
            if (argn > 1)
                sscanf(args[1], "%ld", &xres);
            if (argn > 2)
                sscanf(args[2], "%ld", &yres);
            if (argn > 3)
                sscanf(args[3], "%d", &dang);
            if (argn > 4)
                sscanf(args[4], "%g", &aspect);
            zres = xres;
            break;

        case 11:
            fprintf(stderr, "main - unknown command\n");
            break;

        case -1:
            fprintf(stderr, "main - ambiguous command\n");
            break;
        }

        if (field[0]) {
            MAKERXSTRING(result.shvvalue, field, strlen(field));
            RexxVariablePool(&result);
        }
        returnstring->strptr = NULL;
        returnstring->strlength = 0;
        if (args0)
            free(args0);
    }
    return 0;
}

int main(int argc, char *argv[]) {
    int rc;
    short returnCode;
    RXSTRING Result;

    if (argc > 1) {
        made = 0;
        reset();
        rc = RexxRegisterSubcomExe(hostname, (PFN) planetHandler, NULL);

        Result.strlength = 200;
        Result.strptr = malloc(200);

        rc = RexxStart(0, NULL, argv[1], 0, hostname, RXCOMMAND, NULL,
                &returnCode, &Result);
        if (rc < 0)
            rc = -rc;
    }
    return 0;
}
