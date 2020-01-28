#ifndef POLYLIB_H
#define POLYLIB_H

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/**
 * polylib.h
 *
 * Types and routines for manipulating polynomials
 */

/**
 * Coefficients p[0] to p[N] for a degree-N polynomial
 *   p[0] + p[1]*x + ... + p[N]*x^n
 */
typedef struct Polynomial {
    int N;
    double *p;
} polynomial_t;

/**
 * Pretty-print a polynomial
 */
void poly_print(polynomial_t poly) {
    if (poly.N == 0) {
        printf("%.2e\n", poly.p[0]);
    } else {
        printf("%.2e + ", poly.p[0]);
    }
    for (int i = 1; i < poly.N; i++) {
        printf("%.2e x^%d + ", poly.p[i], i);
    }
    if (poly.N > 0) {
        printf("%.2e x^%d\n", poly.p[poly.N], poly.N);
    }
}

/**
 * Create a polynomial of degree N
 */
polynomial_t poly_new(int N) {
    double *p = malloc((N + 1) * sizeof(double));
    for (int i = 0; i <= N; i++) {
        p[i] = 0.0;
    }
    polynomial_t poly = { N, p };
    return poly;
}

/**
 * Free a polynomial
 */
void poly_free(polynomial_t poly) {
    free(poly.p);
}

/**
 * Evaluate a polynomial
 */
double poly_eval(polynomial_t poly, double x) {
    double val = 0;
    for (int i = 0; i <= poly.N; i++) {
        val += poly.p[i] * pow(x, (double) i);
    }
    return val;
}

/**
 * Multiply a polynomial by its argument
 */
polynomial_t poly_mul_x(polynomial_t poly) {
    polynomial_t poly2 = poly_new(poly.N + 1);
    for (int i = 0; i <= poly.N; i++) {
        poly2.p[i+1] = poly.p[i];
    }
    return poly2;
}

/**
 * Multiply a polynomial by a constant
 */
polynomial_t poly_mul_a(polynomial_t poly, double a) {
    polynomial_t poly2 = poly_new(poly.N);
    for (int i = 0; i <= poly.N; i++) {
        poly2.p[i] = a * poly.p[i];
    }
    return poly2;
}

/**
 * Add two polynomials
 */
polynomial_t poly_add(polynomial_t poly1, polynomial_t poly2) {
    int N_max = poly1.N > poly2.N ? poly1.N : poly2.N;
    polynomial_t poly = poly_new(N_max);
    for (int i = 0; i <= N_max; i++) {
        if (i <= poly1.N && i <= poly2.N) {
            poly.p[i] = poly1.p[i] + poly2.p[i];
        } else if (i <= poly1.N) {
            poly.p[i] = poly1.p[i];
        } else {
            poly.p[i] = poly2.p[i];
        }
    }
    return poly;
}

/**
 * Differentiate a polynomial
 */
polynomial_t poly_diff(polynomial_t poly) {
    polynomial_t poly2 = poly_new(poly.N - 1);
    for (int i = 1; i <= poly.N; i++) {
        poly2.p[i-1] = (double) i * poly.p[i];
    }
    return poly2;
}

/**
 * Create a Legendre polynomial of degree N
 */
polynomial_t legendre_polynomial(int N) {
    polynomial_t poly;
    if (N == 0) {
        poly = poly_new(0);
        poly.p[0] = 1.0;
    } else if (N == 1) {
        poly = poly_new(1);
        poly.p[0] = 0.0;
        poly.p[1] = 1.0;
    } else {
        double n = (double) N;
        polynomial_t p_nminus1 = legendre_polynomial(N - 1);
        polynomial_t p_nminus2 = legendre_polynomial(N - 2);
        polynomial_t x_p_nminus1 = poly_mul_x(p_nminus1);
        polynomial_t c1_x_p_nminus1 = 
            poly_mul_a(x_p_nminus1, (2.0*n - 1.0)/n);
        polynomial_t c2_p_nminus2 = poly_mul_a(p_nminus2, -(n - 1.0)/n);
        poly = poly_add(c1_x_p_nminus1, c2_p_nminus2);
        poly_free(p_nminus1);
        poly_free(p_nminus2);
        poly_free(x_p_nminus1);
        poly_free(c1_x_p_nminus1);
        poly_free(c2_p_nminus2);
    }
    return poly;
}

/**
 * Create a normalized Legendre polynomial of degree N
 */
polynomial_t normalized_legendre_polynomial(int N) {
    polynomial_t leg = legendre_polynomial(N);
    double n = (double) N;
    double norm_factor = sqrt((2.0*n + 1.0)/2.0);
    polynomial_t leg_norm = poly_mul_a(leg, norm_factor);
    poly_free(leg);
    return leg_norm;
}

/**
 * Create a polynomial with roots that define LGL points
 */
polynomial_t lgl_polynomial(int Np) {
    int N = Np - 1;
    polynomial_t p = normalized_legendre_polynomial(N);
    polynomial_t dp = poly_diff(p);
    polynomial_t mdp = poly_mul_a(dp, -1.0);
    polynomial_t xmdp = poly_mul_x(mdp);
    polynomial_t x2mdp = poly_mul_x(mdp);
    polynomial_t lgl_poly = poly_add(dp, x2mdp);
    poly_free(p);
    poly_free(dp);
    poly_free(mdp);
    poly_free(xmdp);
    poly_free(x2mdp);
    return lgl_poly;
}

/**
 * Binary search for polynomial root
 */
double binary_root_search(polynomial_t poly, double left_bound,
        double right_bound, double tol) {

    // Check that bounds surround a root and return NaN if not
    double left_val = poly_eval(poly, left_bound);
    double right_val = poly_eval(poly, right_bound);
    if (left_val * right_val > 0) {
        return NAN;
    }

    // Do binary search for root
    double guess;
    double err;
    do {
        guess = 0.5 * (left_bound + right_bound);
        err = poly_eval(poly, guess);
        if (left_val * err < 0) {
            right_bound = guess;
        } else if (right_val * err < 0) {
            left_bound = guess;
        }
    } while (fabs(err) > tol);

    return guess;
}

/**
 * Calculate LGL quadrature points
 */
double *lgl_quadrature_points(int Np) {

    // Allocate buffer
    double *buf = malloc(Np * sizeof(double));

    // Calculate order N from number of points Np
    int N = Np - 1;

    // Quickly return for some simple cases
    if (N == 0) {
        return buf;
    }
    if (N == 1) {
        buf[0] = -1.0;
        buf[1] = 1.0;
        return buf;
    }
    if (N == 2) {
        buf[0] = -1.0;
        buf[1] = 0.0;
        buf[2] = 1.0;
        return buf;
    }

    // Calculate polynomial with roots that define LGL points
    polynomial_t lgl_poly = lgl_polynomial(Np);

    // Calculate bounds on size of boundary interval
    // See https://arxiv.org/pdf/1311.0028.pdf
    double n = (double) N;
    double min_dr0 = 4.0 / (n*n);
    double max_dr0 = 9.0 * M_PI * M_PI / 
        (2.0 * (2.0*n + 1.0) * (2.0*n + 1.0));

    // Store endpoint nodes
    double node = -1.0;
    int inode = 0;
    buf[inode] = node;
    buf[N - inode] = -node;
    inode++;

    // Search for next node
    double node_min = node + min_dr0;
    double node_max = node + max_dr0;
    double tol = sqrt(DBL_EPSILON);
    node = binary_root_search(lgl_poly, node_min, node_max, tol);
    buf[inode] = node;
    buf[N - inode] = -node;
    inode++;

    // Calculate size of boundary interval
    double dr0 = buf[1] - buf[0];

    // Search for remaining nodes
    while (inode < Np/2) {
        node_min = node + dr0;
        node_max = node_min + dr0;
        while (poly_eval(lgl_poly, node_min) * 
                poly_eval(lgl_poly, node_max) > 0) {
            node_max += dr0;
        }
        node = binary_root_search(lgl_poly, node_min, node_max, tol);
        buf[inode] = node;
        buf[N - inode] = -node;
        inode++;
    }

    // Add center point for even N
    if (N % 2 == 0) {
        buf[N/2] = 0.0;
    }

    poly_free(lgl_poly);
    return buf;
}

#endif
