/* ------------------------------------------------------------
   dVdr_from_rturn.c

   Same mechanism as potscan_cg.c, except:
   - potscan_cg.c scans R at fixed cos(gamma)
   - this code reads R_turn(gamma) from a file
   - evaluates V and derivatives at each given R_turn and angle

   Input file columns:
       gamma_deg   gamma_rad   R_turn_bohr   R_turn_ang

   Output columns:
       gamma_deg gamma_rad cos_gamma R_turn_ang V dVdR dVdr dVdcg abs_dVdr

   Column dVdr = vwrtd = dV/dd, where d is Li2 bond length.
   ------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LCONV .52918
#define ECONV 1.5365e-4
#define DERIVFAC 8.13085e-5
#define EFAC 33.723373
#define VMAX 100000.
#define XINC .05
#define YINC .05
#define XMAX 5.
#define YMAX 5.

#include "potential.h"

double mu = 1.0;
double vr = 1.0;

int main(int argc, char **argv)
{
    FILE *fp;
    char line[512];

    int i;
    double d;
    double v;
    double cg;
    double R;

    double gamma_deg;
    double gamma_rad;
    double R_turn_bohr;
    double R_turn_ang;
    double abs_dVdr;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s rturn_file.txt\n", argv[0]);
        fprintf(stderr, "Input columns expected:\n");
        fprintf(stderr, "gamma_deg gamma_rad R_turn_bohr R_turn_ang\n");
        return 1;
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open input file: %s\n", argv[1]);
        return 1;
    }

    /*
       EXACTLY as in potscan_cg.c
    */
    i = 1;
    d = 3.108;
    Vsetbm();

    printf("# gamma_deg gamma_rad cos_gamma R_turn_ang V dVdR dVdr dVdcg abs_dVdr\n");

    while (fgets(line, sizeof(line), fp) != NULL) {

        if (line[0] == '#') continue;
        if (strlen(line) < 2) continue;
        if (strstr(line, "gamma_deg") != NULL) continue;

        if (sscanf(line, "%lf %lf %lf %lf",
                   &gamma_deg, &gamma_rad, &R_turn_bohr, &R_turn_ang) != 4) {
            continue;
        }

        /*
           Convert angle to cos(gamma), same quantity that potscan_cg.c accepts
           from command line.
        */
        cg = cos(gamma_rad);

        /*
           Avoid tiny numerical roundoff at 90 degrees, e.g. -3e-9 instead of 0.
           This is only a cleanup step.
        */
        if (fabs(cg) < 1.0e-12) {
            cg = 0.0;
        }

        R = R_turn_ang;

        /*
           EXACT same potential call as potscan_cg.c:
               v = V(r, d, cg, i)

           Here r is replaced by the specific R_turn value.
        */
        v = V(R, d, cg, i);

        /*
           EXACT same V scaling as potscan_cg.c.
           Important: potscan only scales V, not vwrtr/vwrtd/vwrtcg.
        */
        v = v * EFAC;

        /*
           EXACT same V cap as potscan_cg.c.
           This only affects the printed potential column, not derivatives.
        */
        if (v > 1000) {
            v = 1000;
        }

        abs_dVdr = fabs(vwrtd);

        printf("%10.4f %16.8f %18.12f %18.8f %18.10e %18.10e %18.10e %18.10e %18.10e\n",
               gamma_deg,
               gamma_rad,
               cg,
               R,
               v,
               vwrtr,
               vwrtd,
               vwrtcg,
               abs_dVdr);
    }

    fclose(fp);
    return 0;
}
