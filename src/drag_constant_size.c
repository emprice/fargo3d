#ifdef CONSTANTDUSTSIZE

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "types_def.h"
#include "define.h"
#include "mpi_dummy.h"
#include "structs.h"
#include "param.h"
#include "fondam.h"
#include "global_ex.h"
#include "prototypes.h"

#ifdef X
#define USEX            (1)
#else
#define USEX            (0)
#endif

#ifdef Y
#define USEY            (1)
#else
#define USEY            (0)
#endif

#ifdef Z
#define USEZ            (1)
#else
#define USEZ            (0)
#endif

#define NDIM            (USEX + USEY + USEZ)

#define MEANMOLWGHT     (2.4)
#define HYDROGENMASS    (1.674e-24 / MSTAR_CGS)
#define CROSSSECTION    (2.367e-15 / (R0_CGS * R0_CGS))

static real pow2(real x) { return x * x; }
static real pow3(real x) { return x * x * x; }

typedef struct
{
    real s;
    real phi;
    real us;
    real uphi;
}
cyl_vector_t;

typedef struct
{
    cyl_vector_t com;
    cyl_vector_t disp;
}
solution_t;

typedef struct
{
    real a;
    real m;
    real rho;
    real cs;
    real nu;
    real lambda_gas;
    const real (*vel_gas)[NDIM];
    const real (*vel_dust)[NDIM];
}
params_t;

static int UpdateSingleRelativeVelocity(real t, params_t *p, real (*dvel)[NDIM]);
static int UpdateRelativeVelocities(real dt, const real (*rho)[NFLUIDS][NDIM],
    const real (*cs)[NDIM], real (*vel)[NFLUIDS][NDIM][NDIM]);

int UpdateSingleRelativeVelocity(real t, params_t *p, real (*dvel)[NDIM])
{
    size_t d;
    real b, fac;

    for (d = 0; d < NDIM; ++d)
    {
        /* save the initial differential velocity */
        (*dvel)[d] = (*(p->vel_dust))[d] - (*(p->vel_gas))[d];
    }

    if (p->lambda_gas > 4 * p->a / 9)
    {
        /* epstein drag regime */
        b = (4 * M_PI * pow2(p->a) * p->rho / (3 * p->m)) *
            sqrt(8 * pow2(p->cs) / (M_PI * GAMMA));
        fac = exp(-b * t);

        /* update velocities */
        for (d = 0; d < NDIM; ++d) (*dvel)[d] *= fac;
    }
    else
    {
        real tf1, tf2;
        real Re, dvmag = 0;

        for (d = 0; d < NDIM; ++d) dvmag += pow2((*dvel)[d]);
        dvmag = sqrt(dvmag);
        Re = 2 * p->a * dvmag / p->nu;

        /* compute the time of first regime transition */
        b = 11 * M_PI * pow2(p->a) * p->rho / (50 * p->m);
        tf1 = (Re - 800) / (800 * b * dvmag);

        if (tf1 > 0)
        {
            /* when tf1 > 0, we know some nonzero time is spent in the
             * large reynolds number regime */
            real t1 = fmin(t, tf1);

            /* large reynolds number regime */
            fac = 1. / (1 + b * t1 * dvmag);

            /* update velocities */
            for (d = 0; d < NDIM; ++d) (*dvel)[d] *= fac;

            /* update time */
            t -= t1;

            /* update magnitude at end of large reynolds number regime */
            Re *= fac;
            dvmag *= fac;
        }

        /* compute the time of second regime transition */
        b = 6 * pow(2, 0.4) * M_PI * pow(p->a, 1.4) * pow(p->nu, 0.6) * p->rho / p->m;
        tf2 = 5 * (pow(Re, 0.4) - 1) / (2 * b * pow(dvmag, 0.4));

        if ((t > 0) && (tf2 > 0))
        {
            /* when tf2 > 0, we know some nonzero time is spent in the
             * intermediate reynolds number regime */
            real t2 = fmin(t, tf2);

            /* intermediate reynolds number regime */
            fac = pow(5. / (5. + 2 * b * t2 * pow(dvmag, 0.4)), 2.5);

            /* update velocities */
            for (d = 0; d < NDIM; ++d) (*dvel)[d] *= fac;

            /* update time */
            t -= t2;

            /* update magnitude at end of intermediate reynolds number regime */
            Re *= fac;
            dvmag *= fac;
        }

        if (t > 0)
        {
            /* any remaining time is spent in the small reynolds number regime */
            b = 6 * M_PI * p->a * p->nu * p->rho / p->m;

            /* small reynolds number regime */
            fac = exp(-b * t);

            /* update velocities */
            for (d = 0; d < NDIM; ++d) (*dvel)[d] *= fac;
        }
    }

    return 0;
}

int UpdateRelativeVelocities(real dt, const real (*rho)[NFLUIDS][NDIM],
    const real (*cs)[NDIM], real (*vel)[NFLUIDS][NDIM][NDIM])
{
#ifndef CYLINDRICAL
#error This code assumes cylindrical coordinates
#endif

    int f, d;
    params_t p;
    real mom, rho_gas, rho_dust, rho_tot;
    real dvel_scratch[NDIM], mu;

    /* update small particle velocities with current gas velocity */
    for (f = 1; f < NFLUIDS; ++f)
    {
        /* skip the one computed separately */
        if (f == BIGDUST) continue;

        p.a = DustRadius[f];
        p.m = DustMass[f];

        for (d = 0; d < NDIM; ++d)
        {
            p.vel_gas = &((*vel)[0][d]);
            p.vel_dust = &((*vel)[f][d]);

            p.cs = (*cs)[d];
            rho_gas = (*rho)[0][d];
            rho_dust = (*rho)[f][d];
            p.rho = rho_gas + rho_dust;

            real mu = ((5 * MEANMOLWGHT * HYDROGENMASS) /
                (64 * CROSSSECTION)) * sqrt(M_PI / GAMMA) * p.cs;
            p.nu = mu / rho_gas;
            p.lambda_gas = sqrt(0.5 * M_PI * GAMMA) * (mu / (rho_gas * p.cs));

            UpdateSingleRelativeVelocity(dt, &p, &dvel_scratch);

            (*vel)[f][d][d] = dvel_scratch[d] + (*vel)[0][d][d];
        }
    }

    /* only one species has feedback */
    {
        f = BIGDUST;

        p.a = DustRadius[f];
        p.m = DustMass[f];

        for (d = 0; d < NDIM; ++d)
        {
            p.vel_gas = &((*vel)[0][d]);
            p.vel_dust = &((*vel)[f][d]);

            p.cs = (*cs)[d];
            rho_gas = (*rho)[0][d];
            rho_dust = (*rho)[f][d];
            p.rho = rho_gas + rho_dust;

            real mu = ((5 * MEANMOLWGHT * HYDROGENMASS) /
                (64 * CROSSSECTION)) * sqrt(M_PI / GAMMA) * p.cs;
            p.nu = mu / rho_gas;
            p.lambda_gas = sqrt(0.5 * M_PI * GAMMA) * (mu / (rho_gas * p.cs));

            mom = rho_gas * (*(p.vel_gas))[d] + rho_dust * (*(p.vel_dust))[d];

            UpdateSingleRelativeVelocity(dt, &p, &dvel_scratch);

            (*vel)[0][d][d] = (mom - rho_dust * dvel_scratch[d]) / p.rho;
            (*vel)[f][d][d] = (mom + rho_gas * dvel_scratch[d]) / p.rho;
        }
    }

    return 0;
}

int SolveDragOde(double dt, unsigned short inplace)
{
#if NFLUIDS > 1
    int i = 0, j = 0, k = 0, f, d;

    real *v_in[NFLUIDS][NDIM];
    real *v_out[NFLUIDS][NDIM];

    const real *Sigma[NFLUIDS];
#ifndef ISOTHERMAL
#error Non-isothermal case not implemented
#endif
#ifdef Z
#error 3-dimensional case not implemented
#endif
    const real *cs = Fluids[0]->Energy->field_cpu;

//<EXTERNAL>
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real* alpha = Alpha;
//<\EXTERNAL>

    for (f = 0; f < NFLUIDS; ++f)
    {
        Sigma[f] = Fluids[f]->Density->field_cpu;
    }

    if (inplace)
    {
        for (f = 0; f < NFLUIDS; ++f)
        {
            d = 0;
#ifdef X
            INPUT(Fluids[f]->Vx_temp);
            OUTPUT(Fluids[f]->Vx_temp);
            v_in[f][d]  = Fluids[f]->Vx_temp->field_cpu;
            v_out[f][d] = Fluids[f]->Vx_temp->field_cpu;
            d++;
#endif
#ifdef Y
            INPUT(Fluids[f]->Vy_temp);
            OUTPUT(Fluids[f]->Vy_temp);
            v_in[f][d]  = Fluids[f]->Vy_temp->field_cpu;
            v_out[f][d] = Fluids[f]->Vy_temp->field_cpu;
            d++;
#endif
        }
    }
    else
    {
        for (f = 0; f < NFLUIDS; ++f)
        {
            d = 0;
#ifdef X
            INPUT(Fluids[f]->Vx);
            OUTPUT(Fluids[f]->Vx_half);
            v_in[f][d]  = Fluids[f]->Vx->field_cpu;
            v_out[f][d] = Fluids[f]->Vx_half->field_cpu;
            d++;
#endif
#ifdef Y
            INPUT(Fluids[f]->Vy);
            OUTPUT(Fluids[f]->Vy_half);
            v_in[f][d]  = Fluids[f]->Vy->field_cpu;
            v_out[f][d] = Fluids[f]->Vy_half->field_cpu;
            d++;
#endif
        }
    }

#ifdef Y
    for (j = 1; j < size_y; ++j)
    {
#endif
#ifdef X
    for (i = XIM; i < size_x; ++i)
    {
#endif
        real soundspeed[NDIM];
        real rho[NFLUIDS][NDIM];
        real vel[NFLUIDS][NDIM][NDIM];

        int lxmyp = lxm + pitch;
        int lxpym = lxp - pitch;

#ifdef X
        {
            real rhom, rhop, csm, csp;
            real hm, hp, Omega;
            size_t d2 = 0;

            csm = cs[lxm]; csp = cs[l];
            soundspeed[d2] = 0.5 * (csm + csp);

            for (f = 0; f < NFLUIDS; ++f)
            {
                Omega = sqrt(G * MSTAR / pow3(Ymed(j)));
                hm = csm / Omega; hp = csp / Omega;

                rhom = Sigma[f][lxm] / (sqrt(2 * M_PI) * hm);
                rhop = Sigma[f][l] / (sqrt(2 * M_PI) * hp);

                rho[f][d2] = 0.5 * (rhom + rhop);

                // x velocity used as-is
                vel[f][d2][0] = v_in[f][0][l];
                // interpolate y-velocity in x
                vel[f][d2][1] = 0.25 * (v_in[f][1][lxm] +
                    v_in[f][1][l] + v_in[f][1][lxmyp] + v_in[f][1][lyp]);
            }

            d2 += 1;
#endif
#ifdef Y
            csm = cs[lym]; csp = cs[l];
            soundspeed[d2] = 0.5 * (csm + csp);

            for (f = 0; f < NFLUIDS; ++f)
            {
                Omega = sqrt(G * MSTAR / pow3(Ymin(j)));
                hm = csm / Omega; hp = csp / Omega;

                rhom = Sigma[f][lym] / (sqrt(2 * M_PI) * hm);
                rhop = Sigma[f][l] / (sqrt(2 * M_PI) * hp);

                rho[f][d2] = 0.5 * (rhom + rhop);

                // interpolate x-velocity in y
                vel[f][d2][0] = 0.25 * (v_in[f][0][lym] +
                    v_in[f][0][l] + v_in[f][0][lxpym] + v_in[f][0][lxp]);
                // y velocity used as-is
                vel[f][d2][1] = v_in[f][1][l];
            }

            d2 += 1;
#endif
        }

        UpdateRelativeVelocities(dt, &rho, &soundspeed, &vel);

        /* copy back to arrays */
        for (f = 0; f < NFLUIDS; ++f)
        {
            for (d = 0; d < NDIM; ++d) v_out[f][d][l] = vel[f][d][d];
        }
#ifdef X
    }
#endif
#ifdef Y
    }
#endif
#endif      /* NFLUIDS > 1 */
    return 0;
}
#endif      /* CONSTANTDUSTSIZE */

/* vim: set ft=c: */
