#ifndef DIAGLIB_H
#define DIAGLIB_H

/**
 * Shallow water model diagnostics
 */

#include "dglib.h"
#include "fieldlib.h"
#include "modellib.h"

/**
 * Potential vorticity diagnostic
 */
typedef struct PVDiagnostic {
    shallowwatermodel_t *swm;
    field2_t u;
    field2_t v;
    field2_t zeta;
    field2_t q;
    field2_t zero;
} pvdiagnostic_t;

/**
 * Create a new PV diagnostic
 */
pvdiagnostic_t pvdiag_new(shallowwatermodel_t *swm) {

    field2_t u = field2_new(swm->h.KX, swm->h.KY, swm->h.dx);
    field2_t v = field2_new(swm->h.KX, swm->h.KY, swm->h.dx);
    field2_t zeta = field2_new(swm->h.KX, swm->h.KY, swm->h.dx);
    field2_t q = field2_new(swm->h.KX, swm->h.KY, swm->h.dx);
    field2_t zero = field2_new(swm->h.KX, swm->h.KY, swm->h.dx);
    field_set_constant(zero, 0.0);

    pvdiagnostic_t pvdiag = { swm, u, v, zeta, q, zero };
    return pvdiag;
}

/**
 * Free memory from a PV diagnostic
 */
void pvdiag_free(pvdiagnostic_t pvdiag) {
    field2_free(pvdiag.u);
    field2_free(pvdiag.v);
    field2_free(pvdiag.zeta);
    field2_free(pvdiag.q);
    field2_free(pvdiag.zero);
}

/**
 * Update PV diagnostic fields
 */
void pvdiag_update(pvdiagnostic_t pvdiag, 
        dgdiff_t diff, dglift_t lift) {

    // Calculate velocity fields
    for (int i = 0; i < pvdiag.swm->h.KX; i++) {
        for (int j = 0; j < pvdiag.swm->h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {
                    pvdiag.u.elements[i][j].mat[ei][ej] =
                        pvdiag.swm->hu.elements[i][j].mat[ei][ej] /
                        pvdiag.swm->h.elements[i][j].mat[ei][ej];
                    pvdiag.v.elements[i][j].mat[ei][ej] =
                        pvdiag.swm->hv.elements[i][j].mat[ei][ej] /
                        pvdiag.swm->h.elements[i][j].mat[ei][ej];
                }
            }
        }
    }

    // Calculate derivatives for vorticity field
    // Note use of PV field as temporary
    // TODO: clean up after refactoring flux computation
    field_set_periodic_boundary_conditions(pvdiag.u);
    field_set_periodic_boundary_conditions(pvdiag.v);
    field_diffx(pvdiag.v, pvdiag.zeta, diff, lift,
            pvdiag.zero, pvdiag.zero, pvdiag.zero, 0.0, 0.0);
    field_diffy(pvdiag.u, pvdiag.q, diff, lift,
            pvdiag.zero, pvdiag.zero, pvdiag.zero, 0.0, 0.0);

    // Finish calculating vorticity and PV
    // Note that for vorticity calculation, q is set to du/dy
    for (int i = 0; i < pvdiag.swm->h.KX; i++) {
        for (int j = 0; j < pvdiag.swm->h.KY; j++) {
            for (int ei = 0; ei < DG_NP; ei++) {
                for (int ej = 0; ej < DG_NP; ej++) {
                    pvdiag.zeta.elements[i][j].mat[ei][ej] -=
                        pvdiag.q.elements[i][j].mat[ei][ej];
                    pvdiag.q.elements[i][j].mat[ei][ej] =
                        (
                         pvdiag.zeta.elements[i][j].mat[ei][ej] +
                         pvdiag.swm->f
                        ) / pvdiag.swm->h.elements[i][j].mat[ei][ej];
                }
            }
        }
    }
}



#endif
