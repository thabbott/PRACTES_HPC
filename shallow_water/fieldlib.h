#ifndef FIELDLIB_H
#define FIELDLIB_H

/**
 * Types and functions for representing fields
 */

#include <stdlib.h>
#include "polylib.h"
#include "lalib.h"
#include "dglib.h"

/**
 * Identifiers for boundaries
 */
typedef enum Direction2 {
    NORTH = 0,
    SOUTH = 1,
    EAST = 2,
    WEST = 3
} direction2_t;

/**
 * Functions on the domain of a 2D field
 */
typedef double (*field2_function_t)(double, double);

/**
 * Two-dimensional field
 */
typedef struct Field2 {
    /**
     * Number of elements in X direction
     */
    int KX;
    /**
     * Number of elements in Y direction
     */
    int KY;
    /**
     * Grid spacing
     */
    double dx;
    /**
     * Interior elements
     */
    dgelement2_t **elements;
    /**
     * Boundary values
     */
    dgboundary2_t **boundaries;
} field2_t;

/**
 * Create a 2D field
 */
field2_t field2_new(int KX, int KY, double dx) {

    // Allocate memory for elements
    dgelement2_t **elements;
    elements = malloc(KX * sizeof(dgelement2_t *));
    for (int i = 0; i < KX; i++) {
        elements[i] = malloc(KY * sizeof(dgelement2_t));
    }

    // Allocate memory for boundaries
    dgboundary2_t **boundaries;
    boundaries = malloc(4 * sizeof(dgboundary2_t *));
    boundaries[NORTH] = malloc(KX * sizeof(dgboundary2_t));
    boundaries[SOUTH] = malloc(KX * sizeof(dgboundary2_t));
    boundaries[EAST] = malloc(KY * sizeof(dgboundary2_t));
    boundaries[WEST] = malloc(KY * sizeof(dgboundary2_t));

    // Create struct and return
    field2_t field = { KX, KY, dx, elements, boundaries };
    return field;
}

/**
 * Free memory for a 2D field
 */
void field2_free(field2_t field) {
    for (int i = 0; i < field.KX; i++) {
        free(field.elements[i]);
    }
    free(field.elements);
    free(field.boundaries[NORTH]);
    free(field.boundaries[SOUTH]);
    free(field.boundaries[EAST]);
    free(field.boundaries[WEST]);
    free(field.boundaries);
}  

/**
 * Fill a pre-allocated matrix from a 2D field. Each entry of the
 * matrix is set to the average value of one of the elements
 * of the field.
 */
void fill_matrix_from_field2(matrix_t mat, field2_t field, 
        dgweights_t weights) {

    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {
            mat.vals[i][j] = dg_element2_average(
                    field.elements[i][j], weights);
        }
    }

}

/**
 * Fill a pre-allocated 2D field from a matrix. Each element
 * is set to a uniform value determined by the entries in
 * the matrix
 */
void fill_field2_from_matrix(field2_t field, matrix_t mat) {

    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {
            field.elements[i][j] = dg_constant_element2(mat.vals[i][j]);
        }
    }

}

/**
 * Return in-element values along an element boundary
 */
dgboundary2_t field_get_element_boundary(
        field2_t field, int i, int j, direction2_t direction) {
    dgboundary2_t boundary;
    switch (direction) {
        case NORTH:
            boundary = dg_ymax_boundary(field.elements[i][j]);
            break;
        case SOUTH:
            boundary = dg_ymin_boundary(field.elements[i][j]);
            break;
        case EAST:
            boundary = dg_xmax_boundary(field.elements[i][j]);
            break;
        case WEST:
            boundary = dg_xmin_boundary(field.elements[i][j]);
            break;
    }
    return boundary;
}

/**
 * Return neighboring values along an element boundary
 */
dgboundary2_t field_get_element_neighbors(
        field2_t field, int i, int j, direction2_t direction) {
    dgboundary2_t boundary;
    switch (direction) {

        case NORTH:
            if (j >= field.KY - 1) {
                boundary = field.boundaries[NORTH][i];
            } else {
                boundary = field_get_element_boundary(
                        field, i, j + 1, SOUTH);
            }
            break;

        case SOUTH:
            if (j <= 0) {
                boundary = field.boundaries[SOUTH][i];
            } else {
                boundary = field_get_element_boundary(
                        field, i, j - 1, NORTH);
            }
            break;

        case EAST:
            if (i >= field.KX - 1) {
                boundary = field.boundaries[EAST][j];
            } else {
                boundary = field_get_element_boundary(
                        field, i + 1, j, WEST);
            }
            break;

        case WEST:
            if (i <= 0) {
                boundary = field.boundaries[WEST][j];
            } else {
                boundary = field_get_element_boundary(
                        field, i - 1, j, EAST);
            }
            break;
    }
    return boundary;
}

/**
 * Set periodic boundary conditions (serial version)
 */
void field_set_periodic_boundary_conditions(field2_t field) {

    int i; int j;

    // North boundary
    j = field.KY - 1;
    for (i = 0; i < field.KX; i++) {
        field.boundaries[NORTH][i] = 
            field_get_element_boundary(field, i, 0, SOUTH);
    }

    // South boundary
    j = 0;
    for (i = 0; i < field.KX; i++) {
        field.boundaries[SOUTH][i] = 
            field_get_element_boundary(field, i, field.KY - 1, NORTH);
    }

    // East boundary
    i = field.KX - 1;
    for (j = 0; j < field.KY; j++) {
        field.boundaries[EAST][j] =
            field_get_element_boundary(field, 0, j, WEST);
    }

    // West boundary
    i = 0;
    for (j = 0; j < field.KY; j++) {
        field.boundaries[WEST][j] =
            field_get_element_boundary(field, field.KX - 1, j, EAST);
    }

}

/**
 * Set one field to the x-derivative of another
 */
void field_diffx(
        field2_t field, field2_t deriv, dgdiff_t diff, dglift_t lift,
        field2_t prognostic, field2_t u, field2_t h, double g,
        double dissipation_factor) {

    dgboundary2_t boundary_low;
    dgboundary2_t boundary_high;
    dgboundary2_t prognostic_low;
    dgboundary2_t prognostic_high;
    dgboundary2_t u_low;
    dgboundary2_t u_high;
    dgboundary2_t h_low;
    dgboundary2_t h_high;
    dgboundary2_t flux_in;
    dgboundary2_t flux_out;
    dgboundary2_t lambda;

    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {

            // Calculate numerical fluxes
            boundary_low = field_get_element_neighbors(
                    field, i, j, WEST);
            boundary_high = field_get_element_boundary(
                    field, i, j, WEST);
            prognostic_low = field_get_element_neighbors(
                    prognostic, i, j, WEST);
            prognostic_high = field_get_element_boundary(
                    prognostic, i, j, WEST);
            u_low = field_get_element_neighbors(
                    u, i, j, WEST);
            u_high = field_get_element_boundary(
                    u, i, j, WEST);
            h_low = field_get_element_neighbors(
                    h, i, j, WEST);
            h_high = field_get_element_boundary(
                    h, i, j, WEST);
            lambda = dg_lax_friedrichs_lambda(
                    u_low, u_high, h_low, h_high, g,
                    dissipation_factor);
            flux_in = dg_lax_friedrichs_flux(
                    boundary_low, boundary_high,
                    prognostic_low, prognostic_high, lambda);
            
            boundary_low = field_get_element_boundary(
                    field, i, j, EAST);
            boundary_high = field_get_element_neighbors(
                    field, i, j, EAST);
            prognostic_low = field_get_element_boundary(
                    prognostic, i, j, EAST);
            prognostic_high = field_get_element_neighbors(
                    prognostic, i, j, EAST);
            u_low = field_get_element_boundary(
                    u, i, j, EAST);
            u_high = field_get_element_neighbors(
                    u, i, j, EAST);
            h_low = field_get_element_boundary(
                    h, i, j, EAST);
            h_high = field_get_element_neighbors(
                    h, i, j, EAST);
            lambda = dg_lax_friedrichs_lambda(
                    u_low, u_high, h_low, h_high, g,
                    dissipation_factor);
            flux_out = dg_lax_friedrichs_flux(
                    boundary_low, boundary_high,
                    prognostic_low, prognostic_high, lambda);

            // Calculate derivatives
            deriv.elements[i][j] = dg_diffx(field.elements[i][j],
                    flux_out, flux_in, diff, lift, 2.0/field.dx);

        }
    }
}

/**
 * Set one field to the y-derivative of another
 */
void field_diffy(
        field2_t field, field2_t deriv, dgdiff_t diff, dglift_t lift,
        field2_t prognostic, field2_t v, field2_t h, double g,
        double dissipation_factor) {

    dgboundary2_t boundary_low;
    dgboundary2_t boundary_high;
    dgboundary2_t prognostic_low;
    dgboundary2_t prognostic_high;
    dgboundary2_t v_low;
    dgboundary2_t v_high;
    dgboundary2_t h_low;
    dgboundary2_t h_high;
    dgboundary2_t flux_in;
    dgboundary2_t flux_out;
    dgboundary2_t lambda;

    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {

            // Calculate numerical fluxes
            boundary_low = field_get_element_neighbors(
                    field, i, j, SOUTH);
            boundary_high = field_get_element_boundary(
                    field, i, j, SOUTH);
            prognostic_low = field_get_element_neighbors(
                    prognostic, i, j, SOUTH);
            prognostic_high = field_get_element_boundary(
                    prognostic, i, j, SOUTH);
            v_low = field_get_element_neighbors(
                    v, i, j, SOUTH);
            v_high = field_get_element_boundary(
                    v, i, j, SOUTH);
            h_low = field_get_element_neighbors(
                    h, i, j, SOUTH);
            h_high = field_get_element_boundary(
                    h, i, j, SOUTH);
            lambda = dg_lax_friedrichs_lambda(
                    v_low, v_high, h_low, h_high, g,
                    dissipation_factor);
            flux_in = dg_lax_friedrichs_flux(
                    boundary_low, boundary_high,
                    prognostic_low, prognostic_high, lambda);
            
            boundary_low = field_get_element_boundary(
                    field, i, j, NORTH);
            boundary_high = field_get_element_neighbors(
                    field, i, j, NORTH);
            prognostic_low = field_get_element_boundary(
                    prognostic, i, j, NORTH);
            prognostic_high = field_get_element_neighbors(
                    prognostic, i, j, NORTH);
            v_low = field_get_element_boundary(
                    v, i, j, NORTH);
            v_high = field_get_element_neighbors(
                    v, i, j, NORTH);
            h_low = field_get_element_boundary(
                    h, i, j, NORTH);
            h_high = field_get_element_neighbors(
                    h, i, j, NORTH);
            lambda = dg_lax_friedrichs_lambda(
                    v_low, v_high, h_low, h_high, g,
                    dissipation_factor);
            flux_out = dg_lax_friedrichs_flux(
                    boundary_low, boundary_high,
                    prognostic_low, prognostic_high, lambda);

            // Calculate derivatives
            deriv.elements[i][j] = dg_diffy(field.elements[i][j],
                    flux_out, flux_in, diff, lift, 2.0/field.dx);

        }
    }
}

/**
 * Apply an in-place filter to a field
 */
void field_filter_inplace(field2_t field, dgfilt_t filter) {
    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {
            field.elements[i][j] = 
                dg_apply_filter(field.elements[i][j], filter);
        }
    }
}

/**
 * Set a field to a constant value
 */
void field_set_constant(field2_t field, double val) {
    for (int i = 0; i < field.KX; i++) {
        for (int j = 0; j < field.KY; j++) {
            field.elements[i][j] = dg_constant_element2(val);
        }
    }
}

/**
 * Set a field from a function
 */
void field_fill_from_function(field2_t field, field2_function_t func) {
    double xlow; double ylow; double x; double y;
    double *r = lgl_quadrature_points(DG_NP);
    for (int i = 0; i < field.KX; i++) {
        xlow = field.dx * (double) i;
        for (int j = 0; j < field.KY; j++) {
            ylow = field.dx * (double) j;
            for (int ei = 0; ei < DG_NP; ei++) {
                x = xlow + 0.5 * (r[ei] - 1.0) * field.dx;
                for (int ej = 0; ej < DG_NP; ej++) {
                    y = ylow + 0.5 * (r[ej] - 1.0) * field.dx;
                    field.elements[i][j].mat[ei][ej] = func(x, y);
                }
            }
        }
    }
    free(r);
}


#endif
