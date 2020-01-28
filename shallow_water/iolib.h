#ifndef IOLIB_H
#define IOLIB_H

/**
 * Routines for I/O
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <png.h>
#include "lalib.h"

/**
 * Convert double to PNG byte
 */
png_byte io_png_byte_from_double(double val, double vmin, double vmax) {
    
    // Rescale input value
    val = (val - vmin)/(vmax - vmin);
    if (val < 0.0) {
        val = 0.0;
    }
    if (val > 1.0) {
        val = 1.0;
    }

    // Convert to png byte
    int ival = (int) (255.0 * val);
    png_byte result = (png_byte) (ival % 256);
    return result;

}

/**
 * Read PNG file as grayscale matrix
 */
int read_png_as_grayscale_matrix(char *filename, matrix_t *mat) {

    // Set up PNG IO and check for errors
    FILE *fp = NULL;
    png_structp png = NULL;
    png_infop info = NULL;
    png_bytep *rows = NULL;
    bool success = true;

    fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Could not open file %s for writing\n",
               filename);
        success = false;
    }
    png = png_create_read_struct(
            PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error creating PNG read struct\n");
        success = false;
    }
    info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error creating PNG info struct\n");
        success = false;
    }
    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error setting jump point\n");
        success = false;
    }

    if (!success) {
        if (fp) fclose(fp);
        if (png || info) png_destroy_read_struct(&png, &info, NULL); 
        return 1;
    } else {
        png_init_io(png, fp);
    }

    // Read header
    png_read_info(png, info);
    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);

    // Strip to 8 bit depth
    if (bit_depth == 16) {
        png_set_strip_16(png);
    }

    // Convert palette to RGB
    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png);
    }

    // Expand gray to 8 bits
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png);
    }

    // Remove alpha channel if present
    if (color_type & PNG_COLOR_MASK_ALPHA) {
        png_set_strip_alpha(png);
    }

    // Convert RGB to grayscale
    if (color_type == PNG_COLOR_TYPE_RGB ||
        color_type == PNG_COLOR_TYPE_RGBA) {
        png_set_rgb_to_gray_fixed(png, 1, -1, -1);
    }

    // Update info
    png_read_update_info(png, info);

    // Allocate row pointers
    rows = malloc(height * sizeof(png_bytep));
    for (int i = 0; i < height; i++) {
        rows[i] = malloc(png_get_rowbytes(png, info));
    }

    // Read image
    png_read_image(png, rows);

    // Clean up
    fclose(fp);
    png_destroy_read_struct(&png, &info, NULL);

    // Create matrix
    *mat = mat_new(height, width);

    // Write values into matrix
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            mat->vals[i][j] = (double) rows[i][j] / 255.0;
        }
    }

    // Free row pointers
    for (int i = 0; i < height; i++) {
        free(rows[i]);
    }
    free(rows);

    return 0;

}

/**
 * Write matrix as grayscale PNG
 */
int write_matrix_as_grayscale_png(char *filename, matrix_t *mat,
        double vmin, double vmax) {

    // Set scale factors if provided values are infinite
    if (isinf(vmin)) {
        vmin = mat_minimum(*mat);
    }
    if (isinf(vmax)) {
        vmax = mat_maximum(*mat);
    }

    // Set up I/O and check for errors
    FILE *fp = NULL;
    png_structp png = NULL;
    png_infop info = NULL;
    png_bytep *rows = NULL;
    bool success = true;

    fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Could not open file %s for writing\n",
               filename);
        success = false;
    }
    png = png_create_write_struct(
            PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error creating PNG write struct\n");
        success = false;
    }
    info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error creating PNG info struct\n");
        success = false;
    }
    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error setting jump point\n");
        success = false;
    }
    if (!success) {
        if (fp) fclose(fp);
        if (png || info) png_destroy_write_struct(&png, &info); 
        return 1;
    } else {
        png_init_io(png, fp);
    }

    // Write header
    int width = mat->N;
    int height = mat->M;
    png_set_IHDR(png, info, width, height, 8,
            PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    // Allocate row pointers
    rows = malloc(height * sizeof(png_bytep));
    for (int i = 0; i < height; i++) {
        rows[i] = malloc(png_get_rowbytes(png, info));
    }

    // Read values from matrix
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            rows[i][j] = io_png_byte_from_double(
                    mat->vals[i][j], vmin, vmax);
        }
    }

    // Write image
    png_write_image(png, rows);
    png_write_end(png, NULL);

    // Free memory
    fclose(fp);
    png_destroy_write_struct(&png, &info);
    for (int i = 0; i < height; i++) {
        free(rows[i]);
    }
    free(rows);

    return 0;

}


#endif
