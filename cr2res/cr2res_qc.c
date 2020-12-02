/*
 * This file is part of the CR2RES Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <string.h>
#include <math.h>

#include <cpl.h>

#include "cr2res_qc.h"
#include "cr2res_trace.h"
#include "cr2res_dfs.h"
#include "cr2res_extract.h"
#include "cr2res_calib.h"
#include "cr2res_detlin.h"
#include "cr2res_wave.h"

#define CR2RES_NONLIN_LEVEL 15000
#define CR2RES_QC_ORDER 4
#define CR2RES_QC_TRACE 1
#define CR2RES_QC_SIZE  100

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_qc  QC related functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The Read Out Noise computation
  @param    ima1        the first input image
  @param    ima2        the second input image
  @param    hsize
  @param    nsamples
  @param    ndit        the NDIT
  @return   the RON or -1 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_dark_qc_ron(
        const cpl_image     *   ima1,
        const cpl_image     *   ima2,
        int                     hsize,
        int                     nsamples,
        int                     ndit)
{
    cpl_image       *   ima ;
    double              norm, ron ;

    /* Test entries */
    if (ima1 == NULL || ima2 == NULL || ndit < 1)   return -1.0 ;

    /* Compute norm */
    norm = 0.5 * ndit ;
    norm = sqrt(norm) ;

    /* Subtraction */
    if ((ima = cpl_image_subtract_create(ima2, ima1)) == NULL) return -1.0 ;

    /* RON measurement */
    cpl_flux_get_noise_window(ima, NULL, hsize, nsamples, &ron, NULL) ;
    cpl_image_delete(ima) ;
    return norm*ron ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin median non linearity 
  @param    coeffs  The detector non linearity coefficients
  @return   The computed median or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_detlin_median(
        const cpl_imagelist     *   coeffs)
{
    double      qc_detlin_median ;
    double      level;
    int         nimgs = 3, width, height;
    hdrl_imagelist * hdrl_coeffs;
    hdrl_image * img;
    hdrl_value value = {CR2RES_NONLIN_LEVEL, 0};
    /* Check Entries */
    if (coeffs == NULL) return -1.0 ;

    /* TODO FIX - > the function returns occasioanlly NaN */
    return 0.0 ;

    /* Initialise */
    qc_detlin_median = -1.0 ;

    // Apply detlin correction on an image with constant value
    hdrl_coeffs = hdrl_imagelist_create((cpl_imagelist*) coeffs, NULL);
    width = hdrl_image_get_size_x(hdrl_imagelist_get(hdrl_coeffs, 0));
    height = hdrl_image_get_size_y(hdrl_imagelist_get(hdrl_coeffs, 0));
    img = hdrl_image_new(width, height); 
    hdrl_image_add_scalar(img, value);
    cr2res_detlin_correct(img, hdrl_coeffs);

    // Then determine the median of that corrected image
    qc_detlin_median = cpl_image_get_median(hdrl_image_get_image(img));

    hdrl_imagelist_delete(hdrl_coeffs);
    hdrl_image_delete(img);

    return qc_detlin_median ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin gain
  @param    coeffs  The detector non linearity coefficients
  @return   The computed gain or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_detlin_gain(
        const cpl_imagelist     *   coeffs)
{
    double      qc_detlin_gain ;

    /* Check Entries */
    if (coeffs == NULL) return -1.0 ;

    /* Initialise */
    qc_detlin_gain = -1.0 ;

    /* TODO */
    
    return qc_detlin_gain ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin min and max level
  @param    ima         input image
  @param    min_level   [out] The computed min level
  @param    max_level   [out] The computed max level
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_detlin_min_max_level(
        const cpl_image     *   ima,
        double              *   min_level,
        double              *   max_level)
{
    double      min_level_loc, max_level_loc ;

    /* Check Entries */
    if (ima == NULL || min_level == NULL || max_level == NULL) return -1 ;

    /* Initialise */
    *min_level = -1.0 ;
    *max_level = -1.0 ;

    /* TODO */
    
    *min_level = min_level_loc ;
    *max_level = max_level_loc ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the central orders positions
  @param    tw          trave Wave table
  @param    order_nb    [out] Array of orders numbers
  @param    order_pos   [out] Array of orders positions
  @param    nbvals      [out] Size of the arrays
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_order_positions(
        const cpl_table *   tw,
        int             **  order_nb,
        double          **  order_pos,
        int             *   nbvals)
{
    const cpl_array *   slit_frac ;
    int             *   order_nb_loc ;
    double          *   order_pos_loc ;
    int                 i, n_full, nrows ;

    /* Check Entries */
    if (tw == NULL || order_nb==NULL || order_pos==NULL || nbvals==NULL) 
        return -1 ;

    /* Initialise */
    *nbvals = 0 ;
    *order_nb = NULL ;
    *order_pos = NULL ;
    nrows = cpl_table_get_nrow(tw) ;
    if (nrows <= 0) return 0 ;
    
    /* Count the number of full slit orders in tw */
    n_full = 0 ;
    for (i=0 ; i<nrows ; i++) {
        slit_frac = cpl_table_get_array(tw, CR2RES_COL_SLIT_FRACTION, i) ;
        if (cr2res_trace_slit_fraction_info(slit_frac, 
                    NULL) == CR2RES_DECKER_NONE) n_full++ ;
    }

    /* Allocate output arrays */
    order_nb_loc = cpl_malloc(n_full * sizeof(int)) ;
    order_pos_loc = cpl_malloc(n_full * sizeof(double)) ;

    /* Fill the arrays */
    n_full = 0 ;
    for (i=0 ; i<nrows ; i++) {
        slit_frac = cpl_table_get_array(tw, CR2RES_COL_SLIT_FRACTION, i) ;
        if (cr2res_trace_slit_fraction_info(slit_frac, 
                    NULL) == CR2RES_DECKER_NONE) {
            order_nb_loc[n_full] = cpl_table_get(tw, CR2RES_COL_ORDER, i, NULL);
            order_pos_loc[n_full] = cr2res_trace_get_trace_ypos(tw, i) ;
            n_full++ ;
        }
    }

    /* Return */
    *nbvals = n_full ;
    *order_nb = order_nb_loc ;
    *order_pos = order_pos_loc ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the total lamp intensity over a limited spectral region
  @param    ima     input image
  @return   The computed intensity or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_lamp_ints(
        const cpl_image     *   ima)
{
    double      qc_lamp_ints ;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_lamp_ints = cpl_image_get_absflux(ima);
    
    return qc_lamp_ints ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean relative background level
  @param    ima     input image
  @return   The computed level or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_mean_level(
        const cpl_image     *   ima)
{
    double      qc_mean_level ;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_mean_level = cpl_image_get_mean(ima);

    /* TODO */
    return qc_mean_level ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the median SNR per pixel
  @param    ima     input image
  @return   The computed median SNR or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_med_snr(
        const cpl_image     *   ima)
{
    double      qc_med_snr ;
    double      sigma;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_med_snr = cpl_image_get_median_dev(ima, &sigma);

    /* TODO */
    
    return qc_med_snr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean and median flux
  @param    ima     input image
  @param    mean_flux   [out] The computed mean flux
  @param    med_flux    [out] The computed median flux
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_mean_med_flux(
        const cpl_image     *   ima,
        double              *   mean_flux,
        double              *   med_flux)
{
    double      mean_flux_loc, med_flux_loc ;

    /* Check Entries */
    if (ima == NULL || mean_flux == NULL || med_flux == NULL) return -1 ;

    /* Initialise */
    *mean_flux = cpl_image_get_mean(ima);
    *med_flux = cpl_image_get_median(ima);
    /* TODO */

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean Y coord of the central order 
  @param    trace   the trace table
  @return   The computed Y coordinate or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_trace_center_y(
        const cpl_table     *   trace)
{
    cpl_vector * vector;
    cpl_array * array;
    cpl_polynomial * poly;
    int * order_idx_values, nb_order_idx_values, central_order_idx, i;
    int * traces, nb_traces;
    double      qc_trace_center_y ;

    /* Check Entries */
    if (trace == NULL) return -1.0 ;

    /* Initialise */
    qc_trace_center_y = 0;
    // Step 1: find central order
    order_idx_values = cr2res_trace_get_order_idx_values((cpl_table*) trace, 
            &nb_order_idx_values);
    array = cpl_array_wrap_int(order_idx_values, nb_order_idx_values);

/* TODO : Is the median really the CENTRAL order ?? */
    central_order_idx = cpl_array_get_median(array);
    cpl_array_unwrap(array);

    // Step 2: Sum all traces together
    traces = cr2res_get_trace_numbers(trace, central_order_idx, &nb_traces);
    for(cpl_size i = 0; i < nb_traces; i++) {
      vector = cr2res_trace_get_ycen(trace, central_order_idx, traces[i], 
                CR2RES_DETECTOR_SIZE);
      qc_trace_center_y += cpl_vector_get_mean(vector);
      cpl_vector_delete(vector);
    }
    
    // Step 3: take the mean
    qc_trace_center_y /= nb_traces;

    cpl_free(order_idx_values);
    cpl_free(traces);
    cpl_polynomial_delete(poly);
    cpl_array_delete(array);

    return qc_trace_center_y ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the number of overexposed pixels in the first raw frame
  @param    ima     the first raw image
  @return   the computed number or -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_nb_overexposed(
        const cpl_image     *   ima)
{
    cpl_mask * mask;
    int     qc_overexposed ;
    int     threshold = CR2RES_DETECTOR_OVEREXP_THRESH;

    /* Check Entries */
    if (ima == NULL) return -1 ;

    /* Initialise */
    qc_overexposed = -1 ;
    
    mask = cpl_mask_threshold_image_create(ima, 0, threshold);
    qc_overexposed = (int) cpl_mask_count(mask);

    cpl_mask_delete(mask);
    return qc_overexposed ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the integrated flux over part of the spectrum
  @param    extracted   Extracted spectrum table
  @return   the computed signal
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_signal(
        const cpl_table     *   extracted)
{
    double          qc_signal ;
    char        *   colname;
    double      *   data;
    cpl_vector  *   vector, *vector2;
    int             nrows, size;

    /* Check Entries */
    if (extracted == NULL) return -1 ;

    /* Initialise */
    qc_signal = -1.0 ;
    colname = cr2res_dfs_SPEC_colname(CR2RES_QC_ORDER, CR2RES_QC_TRACE);
    data = cpl_table_get_data_double((cpl_table*) extracted, colname);
    cpl_free(colname);

    if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
      // QC order or trace not found in the table
      cpl_error_reset();
      return -1;
    }
    
    nrows = cpl_table_get_nrow(extracted);
    size = CR2RES_QC_SIZE;
    if (nrows < size) size = nrows;

    vector = cpl_vector_wrap(nrows, data);
    vector2 = cpl_vector_extract(vector, (nrows - size)/ 2, (nrows + size) 
                / 2 - 1, 1);
    qc_signal = cpl_vector_get_median(vector2);

    cpl_vector_unwrap(vector);
    cpl_vector_delete(vector2);

    return qc_signal ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the SNR of several spectra
  @param    tw          The TW table
  @param    extracted   The Extracted table
  @param    out_order_idx_values    [out] The order values
  @param    out_nb_order_idx_values [out] The number of order values
  @return   An array of size out_nb_order_idx_values with the SNRs values
 */
/*----------------------------------------------------------------------------*/
double * cr2res_qc_snr(
        const cpl_table *   tw,
        const cpl_table *   extracted,
        int             **  out_order_idx_values,
        int             *   out_nb_order_idx_values)
{
    int             *   order_idx_values ;
    int                 nb_order_idx_values ;
    cpl_bivector    *   my_spec,
                    *   my_spec_err ;
    double          *   pmy_spec_err ;
    cpl_vector      *   my_snr_spec ;
    double          *   snrs ;
    int                 i, j ;

    /* Check Entries */
    if (tw==NULL || extracted==NULL || out_order_idx_values==NULL ||
            out_nb_order_idx_values==NULL) return NULL ;

    /* Get the number of orders from TW table */
    order_idx_values = cr2res_trace_get_order_idx_values(tw,
            &nb_order_idx_values) ;

    /* Allocate the output arrays */
    snrs = cpl_malloc(nb_order_idx_values * sizeof(double)) ;

    /* Loop on the orders */
    for (i=0 ; i<nb_order_idx_values ; i++) {
        /* Compute the SNR */
        if (cr2res_extract_EXTRACT1D_get_spectrum(extracted,
                order_idx_values[i], 1, &my_spec, &my_spec_err) == 0) {
            my_snr_spec = cpl_vector_duplicate(cpl_bivector_get_y(my_spec)) ;
            cpl_bivector_delete(my_spec) ;
            /* Clean the error to avoid division by 0.0 */
            pmy_spec_err = cpl_bivector_get_y_data(my_spec_err) ;
            for (j=0 ; j<cpl_bivector_get_size(my_spec_err) ; j++) {
                if (pmy_spec_err[j] == 0) pmy_spec_err[j] = 1.0 ;
            }
            cpl_vector_divide(my_snr_spec, cpl_bivector_get_y(my_spec_err)) ;
            cpl_bivector_delete(my_spec_err) ;
            snrs[i] = cpl_vector_get_median(my_snr_spec) ;
            cpl_vector_delete(my_snr_spec) ; 
        } else {
            snrs[i] = -1.0 ;
        }
    }

    *out_order_idx_values = order_idx_values ; 
    *out_nb_order_idx_values = nb_order_idx_values ;
    return snrs ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the transmission
  @param    extracted   Extracted spectrum table
  @return   the computed signal
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_transmission(
        const cpl_table     *   extracted)
{
    double      qc_transm ;

    /* Check Entries */
    if (extracted == NULL) return -1 ;

    /* Initialise */
    qc_transm = -1.0 ;
    
    /* TODO */
    
    return qc_transm ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the FWHM of the PSF along the slit
  @param    slitfu   Slit func table
  @return   The computed FWHM
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_slit_psf(
        const cpl_table     *   slitfu)
{
    cpl_array  * col_names ;
    cpl_array * fwhm ;
    cpl_vector * x ;
    cpl_vector * y ;
    cpl_vector * save;
    const double     * data ;

    cpl_size i, j ;
    cpl_fit_mode fit_pars ;
    cpl_error_code err;
    int nrow, ncol ;
    double qc_fwhm, x0, sigma, area, offset;

    /* Check Entries */
    if (slitfu == NULL) return -1 ;

    /* Initialise */
    qc_fwhm = -1.0 ;
    fit_pars = CPL_FIT_CENTROID + CPL_FIT_STDEV + CPL_FIT_AREA;
    offset = 0.;

    ncol = cpl_table_get_ncol(slitfu);
    nrow = cpl_table_get_nrow(slitfu);

    fwhm = cpl_array_new(ncol, CPL_TYPE_DOUBLE);
    for (i = 0; i < ncol; i++) cpl_array_set(fwhm, i, -1);
    // x is just the element number
    x = cpl_vector_new(nrow);
    y = cpl_vector_new(nrow);
    for (i = 0; i < nrow; i++){
      cpl_vector_set(x, i, i);
      cpl_vector_set(y, i, 0);
    }
    // Assuming that there are only slitfunc columns in the table
    col_names = cpl_table_get_column_names(slitfu);
    
    for (i = 0; i < ncol; i++)
    {
      // We need to remove "strange" values, i.e. nan and unreasonably large values
      // otherwise the fit will not work
      // Slitfunc should be normalised to the oversampling rate(?)
      data = cpl_table_get_data_double_const(slitfu, cpl_array_get_string(col_names, i));
      for (j = 0; j < nrow; j++) 
      {
        if (isnan(data[j]) | (data[j] > 1)){
          cpl_vector_set(y, j, 0);
        }else{
          cpl_vector_set(y, j, data[j]);}
      }

      err = cpl_vector_fit_gaussian(x, NULL, y, NULL, fit_pars, &x0, &sigma,
                                    &area, &offset, NULL, NULL, NULL);
      qc_fwhm = 2.355 * sigma; // 2.355 = 2 * sqrt(2 * ln(2))
      cpl_array_set(fwhm, i, qc_fwhm);

      if (err != CPL_ERROR_NONE){
        cpl_error_reset();
        cpl_array_set_invalid(fwhm, i);
      }
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG)
    {
      save = cpl_vector_wrap(ncol, cpl_array_get_data_double(fwhm));
      cpl_vector_save(save, "debug_slitfunc_fwhm.fits", CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
      cpl_vector_unwrap(save);
    }

    qc_fwhm = cpl_array_get_median(fwhm);
    
    /* Free Memory */
    cpl_array_delete(col_names);
    cpl_vector_delete(x);
    cpl_vector_delete(y);
    cpl_array_delete(fwhm);

    return qc_fwhm ;
}
/**@}*/

