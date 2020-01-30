#include <stdio.h>
#include <math.h>

#include "iolib.h"
#include "fieldlib.h"
#include "dglib.h"
#include "modellib.h"
#include "diaglib.h"
#include "paralib.h"

int main(int argc, char **argv) {

    // Simulation parameters
    int KX = 384;                       // Grid points per processor (x)
    int KY = 384;                       // Grid points per processor (y)
    int px = 1;                         // Processors (x)
    int py = 1;                         // Processors (x)
    double g = 9.81;                    // Gravity (m/s)
    double f = 1.2e-3;                  // Coriolis parameter (1/s)
    double dissipation_factor = 1e-3;   // Numerical flux dissipation
                                        // factor (0 = centered,
                                        // 1 = upwind)
    char *expid = "loon_384x384";       // Experiment ID
    double h0 = 0.6;                    // Initial height anomaly (m)
    double H = 6.0;                     // Mean height (m)
    double Ld = sqrt(g*H)/f;            // Deformation radius (m)
    double dx = 10.0*Ld/(KX*px);        // Grid spacing (m)
    double CFL = 0.8;                   // CFL number (nondim)
    double dt = CFL * dx / sqrt(g * H); // Time step (s)

    // Time stepping parameters
    int NT = 50;

    // Timer
    double wtime = 0.0;
    double start = 0.0;

    // Initialize parallel computing environment
    para_env_init(&argc, &argv);
    gridtopology_t topo = gridtopology_new(KX, KY, px, py);
    para_printf("Initialized parallel environment\n");

    // Initialize data structures for a shallow water model
    para_printf("Initializing data structures\n");
    shallowwatermodel_t swm = swm_new(KX, KY, dx, g, f, 
            dissipation_factor);
    ssprk3aux_t aux = ssprk3aux_new(KX, KY, dx);
    // rk4aux_t aux = rk4aux_new(KX, KY, dx);
    // rk4coeff_t rk4 = rk4_new();
    dgweights_t weights = dg_weights();
    dgdiff_t diff = dg_diff();
    dglift_t lift = dg_lift();
    dgfilt_t filter = dg_filt(0, 2.0);
    pvdiagnostic_t pvdiag = pvdiag_new(&swm);

    // Set initial conditions and I/O
    para_printf("Setting initial conditions\n");
    field_set_constant(swm.hu, 0.0);
    field_set_constant(swm.hv, 0.0);
    char input[128];
    sprintf(input, "input/%s.png", expid);
    para_printf("Reading height field from %s\n", input);
    matrix_t gmat; matrix_t mat;
    int code = read_png_as_grayscale_matrix(input, &gmat);
    if (code == 0) {
        para_printf("Successful PNG read\n");
    } else {
        para_printf("PNG read failed\n");
        return 1;
    }
    if (gmat.M != KX*px || gmat.N != KY*py) {
        para_printf("Shallow water model (KX*px = %d, KY*py = %d)\n", 
                KX*px, KY*py);
        para_printf("and input PNG (nx = %d, ny = %d)\n", mat.M, mat.N);
        para_printf("have incompatible sizes\n");
        return 2;
    }
    for (int i = 0; i < KX*px; i++) {
        for (int j = 0; j < KY*py; j++) {
            gmat.vals[i][j] = H + h0 * gmat.vals[i][j];
        }
    }
    mat = gridtopo_matrix_subset(gmat, topo);
    fill_field2_from_matrix(swm.h, mat);

    // Append process ID to experiment ID for output
    char outid[128];
    char output[128];
    sprintf(outid, "%s_%03d", expid, para_pid());
    
    // Save initial state
    sprintf(output, "output/%s_%05d.png", outid, 0);
    para_printf("Writing initial condition to %s\n", output);
    fill_matrix_from_field2(mat, swm.h, weights);
    code = write_matrix_as_grayscale_png(
            output, &mat, H, H + h0);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }
    sprintf(output, "output/%s_zeta_%05d.png", outid, 0);
    para_printf("Writing initial vorticity field to %s\n", output);
    pvdiag_update(pvdiag, diff, lift);
    fill_matrix_from_field2(mat, pvdiag.zeta, weights);
    code = write_matrix_as_grayscale_png(
            output, &mat, -INFINITY, INFINITY);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }
    sprintf(output, "output/%s_pv_%05d.png", outid, 0);
    para_printf("Writing initial PV field to %s\n", output);
    fill_matrix_from_field2(mat, pvdiag.q, weights);
    double q0min = mat_minimum(mat);
    double q0max = mat_maximum(mat);
    double q0range = q0max - q0min;
    double stretch = 0.0;
    q0min -= 0.5 * stretch * q0range;
    q0max += 0.5 * stretch * q0range;
    code = write_matrix_as_grayscale_png(
            output, &mat, q0min, q0max);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }

    // Integrate
    para_printf("Starting integration (%d time steps)\n", NT);
    start = paralib_time();
    for (int it = 1; it <= NT; it++) {
        swm_integrate_ssprk3(swm, aux, diff, lift, filter, dt, topo);
    }
    wtime += paralib_elapsed_time_ms(start);
    para_printf("Finished integration\n");

    // Save final state
    int iprint = NT;
    sprintf(output, "output/%s_%05d.png", outid, iprint);
    para_printf("Writing output to %s\n", output);
    fill_matrix_from_field2(mat, swm.h, weights);
    code = write_matrix_as_grayscale_png(
            output, &mat, H, H + h0);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }
    sprintf(output, "output/%s_zeta_%05d.png", outid, iprint);
    para_printf("Writing vorticity field to %s\n", output);
    pvdiag_update(pvdiag, diff, lift);
    fill_matrix_from_field2(mat, pvdiag.zeta, weights);
    code = write_matrix_as_grayscale_png(
            output, &mat, -INFINITY, INFINITY);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }
    sprintf(output, "output/%s_pv_%05d.png", outid, iprint);
    para_printf("Writing PV field to %s\n", output);
    fill_matrix_from_field2(mat, pvdiag.q, weights);
    code = write_matrix_as_grayscale_png(
            output, &mat, q0min, q0max);
    if (code == 0) {
        para_printf("Successful PNG write\n");
    } else {
        para_printf("PNG write failed\n");
        return 3;
    }
    para_printf("Bounds on final values:\n");
    fill_matrix_from_field2(mat, swm.h, weights);
    para_printf("%.2e < h < %.2e\n", 
            mat_minimum(mat), mat_maximum(mat));
    fill_matrix_from_field2(mat, swm.hu, weights);
    para_printf("%.2e < hu < %.2e\n", 
            mat_minimum(mat), mat_maximum(mat));
    fill_matrix_from_field2(mat, swm.hu, weights);
    para_printf("%.2e < hv < %.2e\n", 
            mat_minimum(mat), mat_maximum(mat));

    // Print timing
    para_synchronize();
    para_printf("Time spent on model integration:\n");
    para_synchronize();
    for (int pid = 0; pid < para_nproc(); pid++) {
        if (para_pid() == pid) {
            printf("Process %d: %d milliseconds\n", 
                    para_pid(), (int) wtime);
        }
        para_synchronize();
    }

    // Free memory
    swm_free(swm);
    ssprk3aux_free(aux);
    mat_free(gmat);
    mat_free(mat);

    // Finalize parallel computing environment
    para_env_finalize();

    return 0;

}
