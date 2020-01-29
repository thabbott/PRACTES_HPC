#ifndef PARALIB_H
#define PARALIB_H

/**
 * Types and procedures for parallelization
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdarg.h>

#include "lalib.h"
#include "dglib.h"
#include "fieldlib.h"
#include "mpi.h"

typedef struct GridTopology {
    int pid;
    int KX;
    int KY;
    int px;
    int py;
} gridtopology_t;

/**
 * Error checking
 */
void para_check_err(int ierr, const char *loc) {
    if (ierr != MPI_SUCCESS) {
        printf("%s: fatal MPI error (code %d)\n", loc, ierr);
        exit(ierr);
    }
}

/**
 * Initialize parallel environment
 */
void para_env_init(int *argc, char ***argv) {
    para_check_err(MPI_Init(argc, argv), "para_env_init");
}

/**
 * Finalize parallel environment
 */
void para_env_finalize() {
    para_check_err(MPI_Finalize(), "para_env_finalize");
}

/**
 * Check if processor is root
 */
bool para_is_root() {
    int id;
    para_check_err(
            MPI_Comm_rank(MPI_COMM_WORLD, &id), 
            "para_is_root");
    return (id == 0);
}

/**
 * Get process ID
 */
int para_pid() {
    int id;
    para_check_err(
            MPI_Comm_rank(MPI_COMM_WORLD, &id),
            "para_pid");
    return id;
}

/**
 * Get number of processes
 */
int para_nproc() {
    int nproc;
    para_check_err(
            MPI_Comm_size(MPI_COMM_WORLD, &nproc),
            "para_nproc");
    return nproc;
}

/**
 * Synchronize all processes
 */
void para_synchronize() {
    para_check_err(
            MPI_Barrier(MPI_COMM_WORLD),
            "para_synchronize");
}

/**
 * Print from root processor only
 */
int para_printf(const char *format, ...) {
    va_list arglist;
    int ierr = 0;
    if (para_is_root()) {
        va_start(arglist, format);
        ierr = vprintf(format, arglist);
        va_end(arglist);
    }
    return ierr;
}

/**
 * Get reference time
 */
double paralib_time() {
    return MPI_Wtime();
}

/**
 * Get elapsed time in milliseconds
 */
double paralib_elapsed_time_ms(double ref) {
    return 1e3 * (MPI_Wtime() - ref);
}

/**
 * Create new processor topology
 */
gridtopology_t gridtopology_new(int KX, int KY, int px, int py) {
    int p;
    para_check_err(
            MPI_Comm_size(MPI_COMM_WORLD, &p),
            "gridtopology_new");
    if (p < px*py) {
        para_printf(
            "gridtopology_new: parallel program is oversubscribed.\n");
        para_printf("px = %d, py = %d, p = %d\n", px, py, p);
        exit(1);
    }
    if (p > px*py) {
        para_printf(
            "gridtopology_new: parallel program is undersubscribed.\n");
        para_printf("px = %d, py = %d, p = %d\n", px, py, p);
        exit(1);
    }
    gridtopology_t topo = { para_pid(), KX, KY, px, py };
    return topo;
}

/**
 * Grid topology x-rank
 */
int gridtopo_xrank(gridtopology_t topo) {
    return topo.pid / topo.py;
}

/**
 * Grid topology y-rank
 */
int gridtopo_yrank(gridtopology_t topo) {
    return topo.pid % topo.py;
}

/**
 * PID from ranks
 */
int gridtopo_pid(gridtopology_t topo, int xrank, int yrank) {
    return yrank + xrank*topo.py;
}

/**
 * Grid topology global x-index (low)
 */
int gridtopo_xlow(gridtopology_t topo) {
    return topo.KX * gridtopo_xrank(topo);
}

/**
 * Grid topology global x-index (high + 1)
 */
int gridtopo_xhigh(gridtopology_t topo) {
    return topo.KX * (gridtopo_xrank(topo) + 1);
}

/**
 * Grid topology global y-index (low)
 */
int gridtopo_ylow(gridtopology_t topo) {
    return topo.KY * gridtopo_yrank(topo);
}

/**
 * Grid topology global y-index (high + 1)
 */
int gridtopo_yhigh(gridtopology_t topo) {
    return topo.KY * (gridtopo_yrank(topo) + 1);
}

/**
 * Return neighbor for send
 */
int gridtopo_send_neighbor(gridtopology_t topo, direction2_t direction) {
    int pid;
    int xrank = gridtopo_xrank(topo);
    int yrank = gridtopo_yrank(topo);
    switch (direction) {
        case NORTH:
            yrank = (yrank > 0 ? yrank - 1 : topo.py - 1);
            pid = gridtopo_pid(topo, xrank, yrank);
            break;
        case SOUTH:
            yrank = (yrank < topo.py - 1 ? yrank + 1 : 0);
            pid = gridtopo_pid(topo, xrank, yrank);
            break;
        case EAST:
            xrank = (xrank > 0 ? xrank - 1 : topo.px - 1);
            pid = gridtopo_pid(topo, xrank, yrank);
            break;
        case WEST:
            xrank = (xrank < topo.px - 1 ? xrank + 1 : 0);
            pid = gridtopo_pid(topo, xrank, yrank);
            break;
    }
    return pid;
}

/**
 * Return neighbor for receive
 */
int gridtopo_recv_neighbor(gridtopology_t topo, direction2_t direction) {
    int pid;
    switch (direction) {
        case NORTH:
            pid = gridtopo_send_neighbor(topo, SOUTH);
            break;
        case SOUTH:
            pid = gridtopo_send_neighbor(topo, NORTH);
            break;
        case EAST:
            pid = gridtopo_send_neighbor(topo, WEST);
            break;
        case WEST:
            pid = gridtopo_send_neighbor(topo, EAST);
            break;
    }
    return pid;
}

/**
 * Return local matrix subset for grid topology
 */
matrix_t gridtopo_matrix_subset(matrix_t mat, gridtopology_t topo) {
    matrix_t subset = mat_new(topo.KX, topo.KY);
    int i = 0; int j = 0;
    printf("Process %d owns (%d, %d), (%d, %d)\n",
            para_pid(),
            gridtopo_xlow(topo),
            gridtopo_xhigh(topo),
            gridtopo_ylow(topo),
            gridtopo_yhigh(topo));
    for (int ig = gridtopo_xlow(topo);
         ig < gridtopo_xhigh(topo);
         ig++) {
        j = 0;
        for (int jg = gridtopo_ylow(topo);
             jg < gridtopo_yhigh(topo);
             jg++) {
            subset.vals[i][j] = mat.vals[ig][jg];
            j++;
        }
        i++;
    }
    return subset;
}

/**
 * Fill a buffer from boundary value before exchange
 */
void gridtopo_fill_buffer(field2_t field, direction2_t direction,
        double *buffer) {
    
    dgboundary2_t boundary;
    int ibuf = 0;

    switch (direction) {
        case NORTH:
            for (int i = 0; i < field.KX; i++) {
                boundary = field_get_element_boundary(
                        field, i, 0, SOUTH);
                for (int k = 0; k < DG_NP; k++) {
                    buffer[ibuf] = boundary.mat[k];
                    ibuf++;
                }
            }
            break;

        case SOUTH:
            for (int i = 0; i < field.KX; i++) {
                boundary = field_get_element_boundary(
                        field, i, field.KY - 1, NORTH);
                for (int k = 0; k < DG_NP; k++) {
                    buffer[ibuf] = boundary.mat[k];
                    ibuf++;
                }
            }
            break;
        case EAST:

            for (int i = 0; i < field.KY; i++) {
                boundary = field_get_element_boundary(
                        field, 0, i, WEST);
                for (int k = 0; k < DG_NP; k++) {
                    buffer[ibuf] = boundary.mat[k];
                    ibuf++;
                }
            }
            break;
        case WEST:

            for (int i = 0; i < field.KY; i++) {
                boundary = field_get_element_boundary(
                        field, field.KX - 1, i, EAST);
                for (int k = 0; k < DG_NP; k++) {
                    buffer[ibuf] = boundary.mat[k];
                    ibuf++;
                }
            }
            break;
    }
}

/**
 * Extract boundary values from buffer after exchange
 */
void gridtopo_extract_buffer(field2_t field, direction2_t direction,
        double *buffer) {
    int K = (direction == NORTH || direction == SOUTH ?
            field.KX : field.KY);
    int ibuf = 0;
    for (int i = 0; i < K; i++) {
        for (int k = 0; k < DG_NP; k++) {
            field.boundaries[direction][i].mat[k] = buffer[ibuf];
            ibuf++;
        }
    }
}

/**
 * Exchange boundary values of 2D field
 */
void gridtopo_exchange_boundary(field2_t field, 
        direction2_t direction, gridtopology_t topo) {

    MPI_Request req[2];
    MPI_Status stat[2];

    // Create buffers
    int K = (direction == NORTH || direction == SOUTH ?
            field.KX : field.KY);
    double recv[K*DG_NP];
    double send[K*DG_NP];

    // Fill send buffer
    gridtopo_fill_buffer(field, direction, send);

    // Exchange
    para_check_err(MPI_Isend(
                send,
                K*DG_NP,
                MPI_DOUBLE,
                gridtopo_send_neighbor(topo, direction),
                direction,
                MPI_COMM_WORLD,
                &(req[0])), "gridtopo_exchange_field2");
    para_check_err(MPI_Irecv(
                recv,
                K*DG_NP,
                MPI_DOUBLE,
                gridtopo_recv_neighbor(topo, direction),
                direction,
                MPI_COMM_WORLD,
                &(req[1])), "gridtopo_exchange_field2");
    para_check_err(MPI_Waitall(
                2, 
                req, 
                stat), "gridtopo_exchange_field2");

    // Extract from receive buffer
    gridtopo_extract_buffer(field, direction, recv);
}

/**
 * Set periodic boundary conditions (parallel version)
 */
void field_exchange_periodic_boundary_conditions(field2_t field, 
        gridtopology_t topo) {

    gridtopo_exchange_boundary(field, NORTH, topo);
    gridtopo_exchange_boundary(field, SOUTH, topo);
    gridtopo_exchange_boundary(field, EAST, topo);
    gridtopo_exchange_boundary(field, WEST, topo);
}

#endif
