#ifndef MODELLIB_H
#define MODELLIB_H

/**
 * Types and procedures for integrating a shallow water model
 */

#include "fieldlib.h"
#include "paralib.h"

/**
 * Shallow water model
 */
typedef struct ShallowWaterModel {
    /**
     * Depth
     */
    field2_t h;
    /**
     * U-momentum
     */
    field2_t hu;
    /**
     * V-momentum
     */
    field2_t hv;
    
    /**
     * U-velocity
     */
    field2_t u;
    /**
     * V-velocity
     */
    field2_t v;
    /**
     * Height flux (u)
     */
    field2_t Fhu;
    /**
     * Height flux divergence (u)
     */
    field2_t dxFhu;
    /**
     * Height flux (v)
     */
    field2_t Fhv;
    /**
     * Height flux divergence (v)
     */
    field2_t dyFhv;
    /**
     * U flux (u)
     */
    field2_t Fhuu;
    /**
     * U flux divergence (u)
     */
    field2_t dxFhuu;
    /**
     * U flux (v)
     */
    field2_t Fhuv;
    /**
     * U flux divergence (v)
     */
    field2_t dyFhuv;
    /**
     * V flux (u)
     */
    field2_t Fhvu;
    /**
     * V flux divergence (u)
     */
    field2_t dxFhvu;
    /**
     * V flux (v)
     */
    field2_t Fhvv;
    /**
     * V flux divergence (v)
     */
    field2_t dyFhvv;

    /**
     * Gravity
     */
    double g;
    /**
     * Coriolis parameter
     */
    double f;
    /**
     * Lax-Friedrichs dissipation factor
     */
    double dissipation_factor;

} shallowwatermodel_t;

/**
 * Tendencies and auxiliary storage for strong-stability-preserving RK3
 */
typedef struct SSPRK3Aux {
    field2_t dth;
    field2_t dthu;
    field2_t dthv;
    field2_t h_orig;
    field2_t hu_orig;
    field2_t hv_orig;
} ssprk3aux_t;

/**
 * Tendencies for low-storage RK4
 */
typedef struct RK4Aux {
    field2_t dth;
    field2_t dthu;
    field2_t dthv;
} rk4aux_t;

/**
 * Low-storage RK4 coefficients
 */
typedef struct RK4Coeff {
    double a[5];
    double b[5];
    double c[5];
} rk4coeff_t;

/**
 * Create a new shallow water model
 */
shallowwatermodel_t swm_new(
        int KX, int KY, double dx, double g, double f, 
        double dissipation_factor) {

    field2_t h = field2_new(KX, KY, dx);
    field2_t hu = field2_new(KX, KY, dx);
    field2_t hv = field2_new(KX, KY, dx);
    field2_t u = field2_new(KX, KY, dx);
    field2_t v = field2_new(KX, KY, dx);
    field2_t Fhu = field2_new(KX, KY, dx);
    field2_t dxFhu = field2_new(KX, KY, dx);
    field2_t Fhv = field2_new(KX, KY, dx);
    field2_t dyFhv = field2_new(KX, KY, dx);
    field2_t Fhuu = field2_new(KX, KY, dx);
    field2_t dxFhuu = field2_new(KX, KY, dx);
    field2_t Fhuv = field2_new(KX, KY, dx);
    field2_t dyFhuv = field2_new(KX, KY, dx);
    field2_t Fhvu = field2_new(KX, KY, dx);
    field2_t dxFhvu = field2_new(KX, KY, dx);
    field2_t Fhvv = field2_new(KX, KY, dx);
    field2_t dyFhvv = field2_new(KX, KY, dx);

    shallowwatermodel_t swm = {
        h, hu, hv,
        u, v, 
        Fhu, dxFhu,
        Fhv, dyFhv,
        Fhuu, dxFhuu,
        Fhuv, dyFhuv,
        Fhvu, dxFhvu,
        Fhvv, dyFhvv,
        g, f, dissipation_factor
    };
    return swm;
}

/**
 * Free memory for a shallow water model
 */
void swm_free(shallowwatermodel_t swm) {
    field2_free(swm.h);
    field2_free(swm.hu);
    field2_free(swm.hv);
    field2_free(swm.u);
    field2_free(swm.v);
    field2_free(swm.Fhu);
    field2_free(swm.dxFhu);
    field2_free(swm.Fhv);
    field2_free(swm.dyFhv);
    field2_free(swm.Fhuu);
    field2_free(swm.dxFhuu);
    field2_free(swm.Fhuv);
    field2_free(swm.dyFhuv);
    field2_free(swm.Fhvu);
    field2_free(swm.dxFhvu);
    field2_free(swm.Fhvv);
    field2_free(swm.dyFhvv);
}

/**
 * Create a new set of RK4 tendency fields
 */
rk4aux_t rk4aux_new(int KX, int KY, double dx) {

    field2_t dth = field2_new(KX, KY, dx);
    field2_t dthu = field2_new(KX, KY, dx);
    field2_t dthv = field2_new(KX, KY, dx);
    rk4aux_t aux = {
        dth,
        dthu,
        dthv
    };
    return aux;
}

/**
 * Free memory for RK4 tendency fields
 */
void rk4aux_free(rk4aux_t aux) {
    field2_free(aux.dth);
    field2_free(aux.dthu);
    field2_free(aux.dthv);
}

/**
 * Create a new set of SSP RK3 storage arrays
 */
ssprk3aux_t ssprk3aux_new(int KX, int KY, double dx) {

    field2_t dth = field2_new(KX, KY, dx);
    field2_t dthu = field2_new(KX, KY, dx);
    field2_t dthv = field2_new(KX, KY, dx);
    field2_t h_orig = field2_new(KX, KY, dx);
    field2_t hu_orig = field2_new(KX, KY, dx);
    field2_t hv_orig = field2_new(KX, KY, dx);
    ssprk3aux_t aux = {
        dth,
        dthu,
        dthv,
        h_orig,
        hu_orig,
        hv_orig
    };
    return aux;
}

/**
 * Free memory for SSP RK3 storage arrays
 */
void ssprk3aux_free(ssprk3aux_t aux) {
    field2_free(aux.dth);
    field2_free(aux.dthu);
    field2_free(aux.dthv);
    field2_free(aux.h_orig);
    field2_free(aux.hu_orig);
    field2_free(aux.hv_orig);
}

/**
 * Return a set of low-storage RK4 coefficients
 */
rk4coeff_t rk4_new() {
    return (rk4coeff_t) {
        {
            0.0,
            - 567301805773.0 / 1357537059087.0,
            - 2404267990393.0 / 2016746695238.0,
            - 3550918686646.0 / 2091501179385.0,
            - 1275806237668.0 / 842570457699.0
        }, 
        {
            1432997174477.0 / 9575080441755.0,
            5161836677717.0 / 13612068292357.0,
            1720146321549.0 / 2090206949498.0,
            3134564353537.0 / 4481467310338.0,
            2277821191437.0 / 14882151754819.0
        },
        {
            0.0,
            1432997174477.0 / 9575080441755.0,
            2526269341429.0 / 6820363962896.0,
            2006345519317.0 / 3224310063776.0,
            2802321613138.0 / 2924317926251.0
        }
    };
}

/**
 * Diagnose velocities
 */
void swm_diagnose_uv(shallowwatermodel_t swm) {
    for (int i = 0; i < swm.hu.KX; i++) {
        for (int j = 0; j < swm.hu.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {
                    swm.u.elements[i][j].mat[ei][ej] =
                        swm.hu.elements[i][j].mat[ei][ej] /
                        swm.h.elements[i][j].mat[ei][ej];
                    swm.v.elements[i][j].mat[ei][ej] =
                        swm.hv.elements[i][j].mat[ei][ej] /
                        swm.h.elements[i][j].mat[ei][ej];
                }
            }
        }
    }
}

/**
 * Calculate fluxes
 */
void swm_calculate_fluxes(shallowwatermodel_t swm) {

    for (int i = 0; i < swm.hu.KX; i++) {
        for (int j = 0; j < swm.hu.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {

                    // Height fluxes
                    swm.Fhu.elements[i][j].mat[ei][ej] =
                        swm.hu.elements[i][j].mat[ei][ej];
                    swm.Fhv.elements[i][j].mat[ei][ej] =
                        swm.hv.elements[i][j].mat[ei][ej];

                    // U fluxes
                    swm.Fhuu.elements[i][j].mat[ei][ej] =
                        swm.hu.elements[i][j].mat[ei][ej] *
                        swm.u.elements[i][j].mat[ei][ej] +
                        0.5 * swm.g *
                        swm.h.elements[i][j].mat[ei][ej] *
                        swm.h.elements[i][j].mat[ei][ej];
                    swm.Fhuv.elements[i][j].mat[ei][ej] =
                        swm.hu.elements[i][j].mat[ei][ej] *
                        swm.v.elements[i][j].mat[ei][ej];

                    // V fluxes
                    swm.Fhvu.elements[i][j].mat[ei][ej] =
                        swm.hv.elements[i][j].mat[ei][ej] *
                        swm.u.elements[i][j].mat[ei][ej];
                    swm.Fhvv.elements[i][j].mat[ei][ej] =
                        swm.hv.elements[i][j].mat[ei][ej] *
                        swm.v.elements[i][j].mat[ei][ej] +
                        0.5 * swm.g *
                        swm.h.elements[i][j].mat[ei][ej] *
                        swm.h.elements[i][j].mat[ei][ej];

                }
            }
        }
    }
}

/**
 * Prepare for flux divergence computation by setting boundary conditions
 * (serial version)
 */
void swm_update_boundaries(shallowwatermodel_t swm) {
    field_set_periodic_boundary_conditions(swm.h);
    field_set_periodic_boundary_conditions(swm.hu);
    field_set_periodic_boundary_conditions(swm.hv);
    field_set_periodic_boundary_conditions(swm.u);
    field_set_periodic_boundary_conditions(swm.v);
    field_set_periodic_boundary_conditions(swm.Fhu);
    field_set_periodic_boundary_conditions(swm.Fhv);
    field_set_periodic_boundary_conditions(swm.Fhuu);
    field_set_periodic_boundary_conditions(swm.Fhuv);
    field_set_periodic_boundary_conditions(swm.Fhvu);
    field_set_periodic_boundary_conditions(swm.Fhvv);
}

/**
 * Prepare for flux divergence computation by setting boundary conditions
 * (parallel version)
 */
void swm_exchange_boundaries(shallowwatermodel_t swm,
        gridtopology_t topo) {
    field_exchange_periodic_boundary_conditions(swm.h, topo);
    field_exchange_periodic_boundary_conditions(swm.hu, topo);
    field_exchange_periodic_boundary_conditions(swm.hv, topo);
    field_exchange_periodic_boundary_conditions(swm.u, topo);
    field_exchange_periodic_boundary_conditions(swm.v, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhu, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhv, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhuu, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhuv, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhvu, topo);
    field_exchange_periodic_boundary_conditions(swm.Fhvv, topo);
}

/**
 * Calculate flux divergences
 */
void swm_calculate_flux_divergences(
        shallowwatermodel_t swm, dgdiff_t diff, dglift_t lift) {

    field_diffx(swm.Fhu, swm.dxFhu, diff, lift, 
            swm.h, swm.u, swm.h, swm.g, swm.dissipation_factor);
    field_diffy(swm.Fhv, swm.dyFhv, diff, lift,
            swm.h, swm.v, swm.h, swm.g, swm.dissipation_factor);
    field_diffx(swm.Fhuu, swm.dxFhuu, diff, lift,
            swm.hu, swm.u, swm.h, swm.g, swm.dissipation_factor);
    field_diffy(swm.Fhuv, swm.dyFhuv, diff, lift,
            swm.hu, swm.v, swm.h, swm.g, swm.dissipation_factor);
    field_diffx(swm.Fhvu, swm.dxFhvu, diff, lift,
            swm.hv, swm.u, swm.h, swm.g, swm.dissipation_factor);
    field_diffy(swm.Fhvv, swm.dyFhvv, diff, lift,
            swm.hv, swm.v, swm.h, swm.g, swm.dissipation_factor);

}

/**
 * Apply filter to prognostic fields
 */
void swm_filter(shallowwatermodel_t swm, dgfilt_t filter) {

    field_filter_inplace(swm.h, filter);
    field_filter_inplace(swm.hu, filter);
    field_filter_inplace(swm.hv, filter);
}

/**
 * Integrate one model time step.
 * Uses a strong-stability-preserving third-order 
 * Runge-Kutta integrator
 */
void swm_integrate_ssprk3(
        shallowwatermodel_t swm, ssprk3aux_t aux, 
        dgdiff_t diff, dglift_t lift, dgfilt_t filter, double dt,
        gridtopology_t topo) {

    // Save initial state
    for (int i = 0; i < swm.h.KX; i++) {
        for (int j = 0; j < swm.h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {
                    aux.h_orig.elements[i][j].mat[ei][ej] =
                        swm.h.elements[i][j].mat[ei][ej];
                    aux.hu_orig.elements[i][j].mat[ei][ej] =
                        swm.hu.elements[i][j].mat[ei][ej];
                    aux.hv_orig.elements[i][j].mat[ei][ej] =
                        swm.hv.elements[i][j].mat[ei][ej];
                }
            }
        }
    }
    
    // ===========
    // First stage
    // ===========

    // Calculate flux divergence terms
    swm_diagnose_uv(swm);
    swm_calculate_fluxes(swm);
    swm_exchange_boundaries(swm, topo);
    // swm_update_boundaries(swm);
    swm_calculate_flux_divergences(swm, diff, lift);

    // Update tendencies and fields
    for (int i = 0; i < swm.h.KX; i++) {
        for (int j = 0; j < swm.h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {

                    // Height tendency
                    aux.dth.elements[i][j].mat[ei][ej] =
                         -swm.dxFhu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhv.elements[i][j].mat[ei][ej];

                    // U tendency
                    aux.dthu.elements[i][j].mat[ei][ej] =
                         -swm.dxFhuu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhuv.elements[i][j].mat[ei][ej] +
                         swm.f * 
                         swm.hv.elements[i][j].mat[ei][ej];
                    
                    // V tendency
                    aux.dthv.elements[i][j].mat[ei][ej] =
                         -swm.dxFhvu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhvv.elements[i][j].mat[ei][ej] +
                         -swm.f * 
                         swm.hu.elements[i][j].mat[ei][ej];

                    // Height field
                    swm.h.elements[i][j].mat[ei][ej] =
                        aux.h_orig.elements[i][j].mat[ei][ej] +
                        dt *
                        aux.dth.elements[i][j].mat[ei][ej];
                    
                    // U field
                    swm.hu.elements[i][j].mat[ei][ej] =
                        aux.hu_orig.elements[i][j].mat[ei][ej] +
                        dt *
                        aux.dthu.elements[i][j].mat[ei][ej];
                    
                    // V field
                    swm.hv.elements[i][j].mat[ei][ej] =
                        aux.hv_orig.elements[i][j].mat[ei][ej] +
                        dt *
                        aux.dthv.elements[i][j].mat[ei][ej];

                }
            }
        }
    }
    // Filter resulting prognostic fields
    swm_filter(swm, filter);

    // ============
    // Second stage
    // ============
    
    // Calculate flux divergence terms
    swm_diagnose_uv(swm);
    swm_calculate_fluxes(swm);
    swm_exchange_boundaries(swm, topo);
    // swm_update_boundaries(swm);
    swm_calculate_flux_divergences(swm, diff, lift);

    // Update tendencies and fields
    for (int i = 0; i < swm.h.KX; i++) {
        for (int j = 0; j < swm.h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {

                    // Height tendency
                    aux.dth.elements[i][j].mat[ei][ej] =
                         -swm.dxFhu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhv.elements[i][j].mat[ei][ej];

                    // U tendency
                    aux.dthu.elements[i][j].mat[ei][ej] =
                         -swm.dxFhuu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhuv.elements[i][j].mat[ei][ej] +
                         swm.f * 
                         swm.hv.elements[i][j].mat[ei][ej];
                    
                    // V tendency
                    aux.dthv.elements[i][j].mat[ei][ej] =
                         -swm.dxFhvu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhvv.elements[i][j].mat[ei][ej] +
                         -swm.f * 
                         swm.hu.elements[i][j].mat[ei][ej];

                    // Height field
                    swm.h.elements[i][j].mat[ei][ej] =
                        0.25*(
                            3.0*aux.h_orig.elements[i][j].mat[ei][ej] +
                            swm.h.elements[i][j].mat[ei][ej] +
                            dt*aux.dth.elements[i][j].mat[ei][ej]
                        );
                    
                    // U field
                    swm.hu.elements[i][j].mat[ei][ej] =
                        0.25*(
                            3.0*aux.hu_orig.elements[i][j].mat[ei][ej] +
                            swm.hu.elements[i][j].mat[ei][ej] +
                            dt*aux.dthu.elements[i][j].mat[ei][ej]
                        );
                    
                    // V field
                    swm.hv.elements[i][j].mat[ei][ej] =
                        0.25*(
                            3.0*aux.hv_orig.elements[i][j].mat[ei][ej] +
                            swm.hv.elements[i][j].mat[ei][ej] +
                            dt*aux.dthv.elements[i][j].mat[ei][ej]
                        );
                }
            }
        }
    }
    // Filter resulting prognostic fields
    swm_filter(swm, filter);
    
    // ============
    // Third stage
    // ============
    
    // Calculate flux divergence terms
    swm_diagnose_uv(swm);
    swm_calculate_fluxes(swm);
    swm_exchange_boundaries(swm, topo);
    // swm_update_boundaries(swm);
    swm_calculate_flux_divergences(swm, diff, lift);

    // Update tendencies and fields
    for (int i = 0; i < swm.h.KX; i++) {
        for (int j = 0; j < swm.h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {

                    // Height tendency
                    aux.dth.elements[i][j].mat[ei][ej] =
                         -swm.dxFhu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhv.elements[i][j].mat[ei][ej];

                    // U tendency
                    aux.dthu.elements[i][j].mat[ei][ej] =
                         -swm.dxFhuu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhuv.elements[i][j].mat[ei][ej] +
                         swm.f * 
                         swm.hv.elements[i][j].mat[ei][ej];
                    
                    // V tendency
                    aux.dthv.elements[i][j].mat[ei][ej] =
                         -swm.dxFhvu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhvv.elements[i][j].mat[ei][ej] +
                         -swm.f * 
                         swm.hu.elements[i][j].mat[ei][ej];

                    // Height field
                    swm.h.elements[i][j].mat[ei][ej] =
                        (1.0/3.0)*(
                            aux.h_orig.elements[i][j].mat[ei][ej] +
                            2.0*swm.h.elements[i][j].mat[ei][ej] +
                            2.0*dt*aux.dth.elements[i][j].mat[ei][ej]
                        );
                    
                    // U field
                    swm.hu.elements[i][j].mat[ei][ej] =
                        (1.0/3.0)*(
                            aux.hu_orig.elements[i][j].mat[ei][ej] +
                            2.0*swm.hu.elements[i][j].mat[ei][ej] +
                            2.0*dt*aux.dthu.elements[i][j].mat[ei][ej]
                        );
                    
                    // V field
                    swm.hv.elements[i][j].mat[ei][ej] =
                        (1.0/3.0)*(
                            aux.hv_orig.elements[i][j].mat[ei][ej] +
                            2.0*swm.hv.elements[i][j].mat[ei][ej] +
                            2.0*dt*aux.dthv.elements[i][j].mat[ei][ej]
                        );
                }
            }
        }
    }
    // Filter resulting prognostic fields
    swm_filter(swm, filter);
}

/**
 * Integrate one model time step.
 * Uses a low-storage fourth-order Runge-Kutta integrator
 */
void swm_integrate_rk4(
        shallowwatermodel_t swm, rk4aux_t aux, rk4coeff_t rk4, 
        dgdiff_t diff, dglift_t lift, dgfilt_t filter, double dt,
        gridtopology_t topo) {

    for (int irk = 0; irk < 5; irk++) {

        // Calculate flux divergence terms
        swm_diagnose_uv(swm);
        swm_calculate_fluxes(swm);
        swm_exchange_boundaries(swm, topo);
        swm_calculate_flux_divergences(swm, diff, lift);

        // Update tendencies and fields
        for (int i = 0; i < swm.h.KX; i++) {
            for (int j = 0; j < swm.h.KY; j++) {
                for (int ei = 0; ei < DG_NP; ei++) {
                    for (int ej = 0; ej < DG_NP; ej++) {

                        // Height tendency
                        aux.dth.elements[i][j].mat[ei][ej] =
                            rk4.a[irk] * 
                            aux.dth.elements[i][j].mat[ei][ej] +
                            dt * 
                            (
                             -swm.dxFhu.elements[i][j].mat[ei][ej] +
                             -swm.dyFhv.elements[i][j].mat[ei][ej]
                            );

                        // U tendency
                        aux.dthu.elements[i][j].mat[ei][ej] =
                            rk4.a[irk] *
                            aux.dthu.elements[i][j].mat[ei][ej] +
                            dt *
                            (
                             -swm.dxFhuu.elements[i][j].mat[ei][ej] +
                             -swm.dyFhuv.elements[i][j].mat[ei][ej] +
                             swm.f * 
                             swm.hv.elements[i][j].mat[ei][ej]
                            );
                        
                        // V tendency
                        aux.dthv.elements[i][j].mat[ei][ej] =
                            rk4.a[irk] *
                            aux.dthv.elements[i][j].mat[ei][ej] +
                            dt *
                            (
                             -swm.dxFhvu.elements[i][j].mat[ei][ej] +
                             -swm.dyFhvv.elements[i][j].mat[ei][ej] +
                             -swm.f * 
                             swm.hu.elements[i][j].mat[ei][ej]
                            );

                        // Height field
                        swm.h.elements[i][j].mat[ei][ej] =
                            swm.h.elements[i][j].mat[ei][ej] +
                            rk4.b[irk] *
                            aux.dth.elements[i][j].mat[ei][ej];
                        
                        // U field
                        swm.hu.elements[i][j].mat[ei][ej] =
                            swm.hu.elements[i][j].mat[ei][ej] +
                            rk4.b[irk] *
                            aux.dthu.elements[i][j].mat[ei][ej];
                        
                        // V field
                        swm.hv.elements[i][j].mat[ei][ej] =
                            swm.hv.elements[i][j].mat[ei][ej] +
                            rk4.b[irk] *
                            aux.dthv.elements[i][j].mat[ei][ej];

                    }
                }
            }
        }
        // Filter resulting prognostic fields
        swm_filter(swm, filter);
    }
}

/**
 * Integrate one model time step.
 * Uses forward Euler integrator
 */
void swm_integrate_euler(
        shallowwatermodel_t swm,
        dgdiff_t diff, dglift_t lift, dgfilt_t filter, double dt,
        gridtopology_t topo) {

    // Calculate flux divergence terms
    swm_diagnose_uv(swm);
    swm_calculate_fluxes(swm);
    swm_exchange_boundaries(swm, topo);
    swm_calculate_flux_divergences(swm, diff, lift);

    // Update fields
    for (int i = 0; i < swm.h.KX; i++) {
        for (int j = 0; j < swm.h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {

                    // Height
                    swm.h.elements[i][j].mat[ei][ej] +=
                        dt * 
                        (
                         -swm.dxFhu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhv.elements[i][j].mat[ei][ej]
                        );

                    // U
                    swm.hu.elements[i][j].mat[ei][ej] +=
                        dt *
                        (
                         -swm.dxFhuu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhuv.elements[i][j].mat[ei][ej] +
                         -swm.f * 
                         swm.hv.elements[i][j].mat[ei][ej]
                        );
                    
                    // V
                    swm.hv.elements[i][j].mat[ei][ej] +=
                        dt *
                        (
                         -swm.dxFhvu.elements[i][j].mat[ei][ej] +
                         -swm.dyFhvv.elements[i][j].mat[ei][ej] +
                         swm.f * 
                         swm.hu.elements[i][j].mat[ei][ej]
                        );
                }
            }
        }
    }

    // Filter resulting prognostic fields
    swm_filter(swm, filter);
}

#endif
