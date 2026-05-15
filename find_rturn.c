#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EFAC         33.723373   /* trajectory energy units → cm-1  */
#define LCONV        0.52918     /* bohr → Angstrom                 */
#define RMAX         15.0        /* start of inward scan (Angstrom) */
#define COARSE_STEP  0.05        /* coarse scan step (Angstrom)     */
#define TOLERANCE    1e-8        /* bisection convergence (Angstrom)*/

#include "potential.h"

double mu = 1, vr = 1;

int main(int argc, char **argv){

    double cg, E_col, d;
    double r, r_prev, v_curr, v_prev;
    double r_lo, r_hi, r_mid, v_mid;
    double R_turn, V_turn;
    int    i = 0;
    int    iter = 0;
    int    found = 0;

    /* ── command line arguments ──────────────────────────────────────────── */
    if(argc < 3){
        printf("Usage: %s gamma_deg E_col\n", argv[0]);
        printf("  gamma_deg = angle in degrees  (0 = end-on, 90 = side-on)\n");
        printf("  E_col     = collision energy in cm-1\n");
        printf("Example: %s 0 2747\n", argv[0]);
        return 1;
    }

    double gamma_deg = atof(argv[1]);
    cg    = cos(gamma_deg * M_PI / 180.0);
    E_col = atof(argv[2]);
    d     = 3.108;           /* r_eq for Li2(A) in Angstrom */

    /* ── initialise potential ─────────────────────────────────────────────── */
    Vsetbm();

    /* ── coarse scan from RMAX inward to find outer bracket ─────────────── */
    /*                                                                         */
    /*  We scan from large R toward small R looking for the FIRST place       */
    /*  where V rises above E_col. This is the outermost turning point —      */
    /*  the physically correct one for scattering.                            */
    /*                                                                         */
    /*  Bracket condition:                                                     */
    /*    r_hi = larger R  where V < E_col  (particle can be here)            */
    /*    r_lo = smaller R where V >= E_col (classically forbidden)           */

    r_prev = RMAX;
    v_prev = V(r_prev, d, cg, i) * EFAC;

    for(r = RMAX - COARSE_STEP; r >= 0.1; r -= COARSE_STEP){

        v_curr = V(r, d, cg, i) * EFAC;

        /* V was below E_col and just rose above: outer bracket found */
        if(v_prev < E_col && v_curr >= E_col){
            r_hi  = r_prev;   /* large R, V < E_col  */
            r_lo  = r;        /* small R, V >= E_col */
            found = 1;
            break;
        }

        r_prev = r;
        v_prev = v_curr;
    }

    if(!found){
        printf("ERROR: no turning point found for gamma=%.1f E_col=%.2f\n",
               gamma_deg, E_col);
        printf("  V(RMAX=%.1f) = %.4f cm-1\n",
               RMAX, V(RMAX, d, cg, i)*EFAC);
        return 1;
    }

    /* ── bisection between r_lo and r_hi ─────────────────────────────────── */
    /*  r_lo: small R, V >= E_col                                             */
    /*  r_hi: large R, V <  E_col                                             */
    /*  Each step halves the interval — no interpolation, exact V evaluation  */

    while((r_hi - r_lo) > TOLERANCE){
        r_mid = 0.5 * (r_lo + r_hi);
        v_mid = V(r_mid, d, cg, i) * EFAC;

        if(v_mid >= E_col){
            r_lo = r_mid;    /* turning point is further out */
        } else {
            r_hi = r_mid;    /* turning point is further in  */
        }
        iter++;
    }

    R_turn = 0.5 * (r_lo + r_hi);
    V_turn = V(R_turn, d, cg, i) * EFAC;

    printf("gamma       = %.6f degrees\n", gamma_deg);
    printf("cosγ        = %.6f\n",   cg);
    printf("E_col       = %.4f cm-1\n", E_col);
    printf("R_turn      = %.8f Angstrom\n", R_turn);
    printf("R_turn      = %.8f bohr\n",     R_turn / LCONV);
    printf("V(R_turn)   = %.6f cm-1\n",     V_turn);
    printf("iterations  = %d\n", iter);

    return 0;
}
