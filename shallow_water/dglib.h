#ifndef DGLIB_H
#define DGLIB_H

/**
 * Discontinuous Galerkin operators and elements
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "polylib.h"
#include "lalib.h"

/**
 * Number of nodes in each dimension within an element.
 * Set at compile-time to allow automatic stack allocation
 * of operators and elements.
 */
#ifndef DG_N
#define DG_N 1
#endif
#define DG_NP (DG_N + 1)

/**
 * Dimensions
 */
typedef enum DGDimension {
    DGXDIM = 0,
    DGYDIM = 1,
    DGZDIM = 2
} dgdimension_t;

/**
 * Two-dimensional DG element
 */
typedef struct DGElement2 {
    double mat[DG_NP][DG_NP];
} dgelement2_t;

/**
 * Two-dimensional DG element boundary
 */
typedef struct DGBoundary2 {
    double mat[DG_NP];
} dgboundary2_t;

/**
 * Statically-sized square matrix
 */
typedef struct DGMat {
    double mat[DG_NP][DG_NP];
} dgmat_t;

/**
 * Vandermonde matrix
 */
typedef dgmat_t dgvand_t;

/**
 * Averaging weights
 */
typedef dgmat_t dgweights_t;

/**
 * Differentiation operator
 */
typedef dgmat_t dgdiff_t;

/**
 * LIFT operator
 */
typedef struct DGLift {
    double mat[DG_NP][2];
} dglift_t;

/**
 * Filter operator
 */
typedef dgmat_t dgfilt_t;

/**
 * Convert DG Vandermonde matrix to generic matrix
 */
matrix_t vand_to_matrix(dgvand_t vand) {
    matrix_t mat = mat_new(DG_NP, DG_NP);
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            mat.vals[i][j] = vand.mat[i][j];
        }
    }
    return mat;
}

/**
 * Convert DG differentiation operator to generic matrix
 */
matrix_t diff_to_matrix(dgdiff_t diff) {
    matrix_t mat = mat_new(DG_NP, DG_NP);
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            mat.vals[i][j] = diff.mat[i][j];
        }
    }
    return mat;
}

/**
 * Convert DG lift operator to generic matrix
 */
matrix_t lift_to_matrix(dglift_t lift) {
    matrix_t mat = mat_new(DG_NP, 2);
    for (int i = 0; i < DG_NP; i++) {
        mat.vals[i][0] = lift.mat[i][0];
        mat.vals[i][1] = lift.mat[i][1];
    }
    return mat;
}

/**
 * Convert DG filter operator to generic matrix
 */
matrix_t filt_to_matrix(dgfilt_t filt) {
    matrix_t mat = mat_new(DG_NP, DG_NP);
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            mat.vals[i][j] = filt.mat[i][j];
        }
    }
    return mat;
}


/**
 * Modal basis function
 */
polynomial_t dg_psi(int Np) {
    return normalized_legendre_polynomial(Np - 1);
}

/**
 * Create a Vandermonde matrix
 */
dgvand_t dg_vand() {
    dgvand_t vand;
    double *lglpts = lgl_quadrature_points(DG_NP);
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            polynomial_t psi = dg_psi(j + 1);
            vand.mat[i][j] = poly_eval(psi, lglpts[i]);
            poly_free(psi);
        }
    }
    free(lglpts);
    return vand;
}

/**
 * Calculate averaging weights
 */
dgweights_t dg_weights() {
    dgweights_t weights;
    dgvand_t vand = dg_vand();
    matrix_t V = vand_to_matrix(vand);
    matrix_t Vinv = la_inverse(V);
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            weights.mat[i][j] = 0.5*Vinv.vals[0][i]*Vinv.vals[0][j];
        }
    }
    mat_free(V);
    mat_free(Vinv);
    return weights;
}


/**
 * Create a differentiation operator
 */
dgdiff_t dg_diff() {

    dgvand_t vand = dg_vand();
    matrix_t Vr = mat_new(DG_NP, DG_NP);
    matrix_t V = vand_to_matrix(vand);
    double *lglpts = lgl_quadrature_points(DG_NP);

    // Calculate Vr
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            polynomial_t psi = dg_psi(j + 1);
            polynomial_t dpsi = poly_diff(psi);
            Vr.vals[i][j] = poly_eval(dpsi, lglpts[i]);
            poly_free(psi);
            poly_free(dpsi);
        }
    }

    // Calculate differentiation matrix
    matrix_t Vinv = la_inverse(V);
    matrix_t VrVinv = la_product(Vr, Vinv);

    // Copy results to differentiation operator
    dgdiff_t diff;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            diff.mat[i][j] = VrVinv.vals[i][j];
        }
    }

    mat_free(Vr); mat_free(V); free(lglpts);
    mat_free(Vinv); mat_free(VrVinv);
    return diff;
}

/**
 * Create a LIFT operator
 */
dglift_t dg_lift() {

    // Calculate lift matrix
    dgvand_t vand = dg_vand();
    matrix_t V = vand_to_matrix(vand);
    matrix_t VT = la_transpose(V);
    matrix_t VVT = la_product(V, VT);
    matrix_t E = mat_zeros(DG_NP, 2);
    E.vals[0][0] = -1.0;
    E.vals[DG_NP - 1][1] = 1.0;
    matrix_t VVTE = la_product(VVT, E);

    // Copy results to lift operator
    dglift_t lift;
    for (int i = 0; i < DG_NP; i++) {
        lift.mat[i][0] = VVTE.vals[i][0];
        lift.mat[i][1] = VVTE.vals[i][1];
    }

    mat_free(V); mat_free(VT); mat_free(VVT); 
    mat_free(E); mat_free(VVTE);
    return lift;
}

/**
 * Calculate diagonal coefficient for DG filter
 */
double df_filter_coeff(double eta, int cutoff, double skew) {
    double etac = (double) cutoff / (double) DG_N;
    if (eta < etac) {
        return 1.0;
    } else {
        double alpha = -log(DBL_EPSILON);
        return exp(-alpha * pow((eta - etac)/(1.0 - etac), skew));
    }
}

/**
 * Create a filter operator
 * Larger cutoff or skew decreases filter strength
 */
dgfilt_t dg_filt(int cutoff, double skew) {

    // Calculate filter matrix
    matrix_t diag = mat_zeros(DG_NP, DG_NP);
    for (int i = 0; i < DG_NP; i++) {
        diag.vals[i][i] = df_filter_coeff(
                (double) i / (double) DG_N, cutoff, skew);
    }
    dgvand_t vand = dg_vand();
    matrix_t V = vand_to_matrix(vand);
    matrix_t Vinv = la_inverse(V);
    matrix_t Vdiag = la_product(V, diag);
    matrix_t VdiagVinv = la_product(Vdiag, Vinv);

    // Copy results to filter operator
    dgfilt_t filter;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            filter.mat[i][j] = VdiagVinv.vals[i][j];
        }
    }

    mat_free(V); mat_free(Vinv); 
    mat_free(Vdiag); mat_free(VdiagVinv);
    return filter;
}

/**
 * Return interior x-derivative term at a single node of an element
 */
double dg_diffx_interior_term(
        dgelement2_t element, int i, int j, dgdiff_t diff) {
    double result = 0.0;
    for (int k = 0; k < DG_NP; k++) {
        result += diff.mat[j][k] * element.mat[k][j];
    }
    return result;
}

/**
 * Return interior y-derivative term at a single node of an element
 */
double dg_diffy_interior_term(
        dgelement2_t element, int i, int j, dgdiff_t diff) {
    double result = 0.0;
    for (int k = 0; k < DG_NP; k++) {
        result += diff.mat[i][k] * element.mat[i][k];
    }
    return result;
}

/**
 * Return boundary x-derivative term at a single node of an element
 */
double dg_diffx_boundary_term(
        dgelement2_t element, int i, int j,
        dgboundary2_t flux_in, dgboundary2_t flux_out,
        dglift_t lift) {
    double result = 
        -lift.mat[j][0] * (element.mat[0][j] - flux_in.mat[j]);
    result += 
        -lift.mat[j][1] * (element.mat[DG_NP-1][j] - flux_out.mat[j]);
    return result;
}

/**
 * Return boundary y-derivative term at a single node of an element
 */
double dg_diffy_boundary_term(
        dgelement2_t element, int i, int j,
        dgboundary2_t flux_in, dgboundary2_t flux_out,
        dglift_t lift) {
    double result = 
        -lift.mat[i][0] * (element.mat[i][0] - flux_in.mat[i]);
    result +=
        -lift.mat[i][1] * (element.mat[i][DG_NP-1] - flux_out.mat[i]);
    return result;
}

/**
 * Return the x-derivative of an element
 */
dgelement2_t dg_diffx(dgelement2_t element,
        dgboundary2_t flux_out, dgboundary2_t flux_in,
        dgdiff_t diff, dglift_t lift, double metric) {
    
    dgelement2_t deriv;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            deriv.mat[i][j] = 
                dg_diffx_interior_term(element, i, j, diff);
            deriv.mat[i][j] +=
                dg_diffx_boundary_term(element, i, j,
                        flux_in, flux_out, lift);
            deriv.mat[i][j] *= metric;
        }
    }
    return deriv;
}

/**
 * Return the y-derivative of an element
 */
dgelement2_t dg_diffy(dgelement2_t element,
        dgboundary2_t flux_out, dgboundary2_t flux_in,
        dgdiff_t diff, dglift_t lift, double metric) {
    
    dgelement2_t deriv;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            deriv.mat[i][j] = 
                dg_diffy_interior_term(element, i, j, diff);
            deriv.mat[i][j] +=
                dg_diffy_boundary_term(element, i, j,
                        flux_in, flux_out, lift);
            deriv.mat[i][j] *= metric;
        }
    }
    return deriv;
}

/**
 * Return the element boundary with the smallest x-value
 */
dgboundary2_t dg_xmin_boundary(dgelement2_t element) {
    dgboundary2_t boundary;
    for (int i = 0; i < DG_NP; i++) {
        boundary.mat[i] = element.mat[0][i];
    }
    return boundary;
}

/**
 * Return the element boundary with the largest x-value
 */
dgboundary2_t dg_xmax_boundary(dgelement2_t element) {
    dgboundary2_t boundary;
    for (int i = 0; i < DG_NP; i++) {
        boundary.mat[i] = element.mat[DG_NP-1][i];
    }
    return boundary;
}

/**
 * Return the element boundary with the smallest y-value
 */
dgboundary2_t dg_ymin_boundary(dgelement2_t element) {
    dgboundary2_t boundary;
    for (int i = 0; i < DG_NP; i++) {
        boundary.mat[i] = element.mat[i][0];
    }
    return boundary;
}

/**
 * Return the element boundary with the largest y-value
 */
dgboundary2_t dg_ymax_boundary(dgelement2_t element) {
    dgboundary2_t boundary;
    for (int i = 0; i < DG_NP; i++) {
        boundary.mat[i] = element.mat[i][DG_NP-1];
    }
    return boundary;
}

/**
 * Calculate the local Lax-Friedrichs lambda parameter
 */
dgboundary2_t dg_lax_friedrichs_lambda(
        dgboundary2_t u_low, dgboundary2_t u_high,
        dgboundary2_t h_low, dgboundary2_t h_high,
        double g, double dissipation_factor) {
    dgboundary2_t lambda;
    double l_low; double l_high;
    for (int i = 0; i < DG_NP; i++) {
        l_low = fabs(u_low.mat[i]) + sqrt(g*h_low.mat[i]);
        l_high = fabs(u_high.mat[i]) + sqrt(g*h_low.mat[i]);
        lambda.mat[i] = dissipation_factor *
            (l_low > l_high ? l_low : l_high);
    }
    return lambda;
}

/**
 * Calculate a local Lax-Friedrichs numerical flux
 */
dgboundary2_t dg_lax_friedrichs_flux(
        dgboundary2_t boundary_low, dgboundary2_t boundary_high,
        dgboundary2_t value_low, dgboundary2_t value_high,
        dgboundary2_t lambda) {
    dgboundary2_t flux;
    for (int i = 0; i < DG_NP; i++) {
        flux.mat[i] =  0.5*(boundary_low.mat[i] + boundary_high.mat[i])
            - 0.5*lambda.mat[i]*(value_high.mat[i] - value_low.mat[i]);
    }
    return flux;
}

/**
 * Calculate a centered numerical flux
 */
dgboundary2_t dg_centered_flux(
        dgboundary2_t boundary_from, dgboundary2_t boundary_to) {
    dgboundary2_t flux;
    for (int i = 0; i < DG_NP; i++) {
        flux.mat[i] =  0.5*(boundary_from.mat[i] + boundary_to.mat[i]);
    }
    return flux;
}

/**
 * Return a filtered element
 */
dgelement2_t dg_apply_filter(dgelement2_t element, dgfilt_t filter) {
    dgelement2_t aux;
    dgelement2_t filtered_element;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            aux.mat[i][j] = 0.0;
            for (int k = 0; k < DG_NP; k++) {
                aux.mat[i][j] +=
                    element.mat[i][k] * filter.mat[j][k];
            }
        }
    }
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            filtered_element.mat[i][j] = 0.0;
            for (int k = 0; k < DG_NP; k++) {
                filtered_element.mat[i][j] +=
                    filter.mat[i][k] * aux.mat[k][j];
            }
        }
    }
    return filtered_element;
}

/**
 * Return a constant element
 */
dgelement2_t dg_constant_element2(double val) {
    dgelement2_t element;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            element.mat[i][j] = val;
        }
    }
    return element;
}

/**
 * Calculate the average value of an element.
 */
double dg_element2_average(dgelement2_t element, dgweights_t weights) {
    double avg = 0.0;
    for (int i = 0; i < DG_NP; i++) {
        for (int j = 0; j < DG_NP; j++) {
            avg += element.mat[i][j] * weights.mat[i][j];
        }
    }
    return avg;
}


#endif
