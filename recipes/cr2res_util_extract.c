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

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_trace.h"
#include "cr2res_slitdec.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_extract"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_table * cr2res_extract_tab_create(
        cpl_vector      **  spectrum,
        cpl_table       *   trace_table) ;
static cpl_table * cr2res_slit_func_tab_create(
        cpl_vector      **  slit_func,
        cpl_table       *   trace_table) ;
static int cr2res_util_extract_create(cpl_plugin *);
static int cr2res_util_extract_exec(cpl_plugin *);
static int cr2res_util_extract_destroy(cpl_plugin *);
static int cr2res_util_extract(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_extract_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"science.fits " CR2RES_COMMAND_LINE "\n"
"trace.fits " CR2RES_COMMAND_LINE "\n"
" The recipe produces the following products:\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_extract 	Optimal Extraction Utility
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list    the plugin list
  @return   0 if everything is ok, 1 otherwise
  @note     Only this function is exported

  Create the recipe instance and make it available to the application using the
  interface.
 */
/*----------------------------------------------------------------------------*/
int cpl_plugin_get_info(cpl_pluginlist * list)
{
    cpl_recipe  *   recipe = cpl_calloc(1, sizeof *recipe );
    cpl_plugin  *   plugin = &recipe->interface;

    if (cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    CR2RES_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "cr2res_util_extract",
                    "Optimal Extraction utility",
                    cr2res_util_extract_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_extract_create,
                    cr2res_util_extract_exec,
                    cr2res_util_extract_destroy)) {
        cpl_msg_error(cpl_func, "Plugin initialization failed");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    if (cpl_pluginlist_append(list, plugin)) {
        cpl_msg_error(cpl_func, "Error adding plugin to list");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Setup the recipe options
  @param    plugin  the plugin
  @return   0 if everything is ok

  Defining the command-line/configuration parameters for the recipe.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;
    cpl_parameter * p;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_util_extract", 10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_util_extract", 64);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_util_extract", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_util_extract", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.sum_only",
            CPL_TYPE_BOOL, "Flag to only sum along detector",
            "cr2res.cr2res_util_extract", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "sum_only");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_extract", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_extract", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_extract", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
 
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_extract(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract_destroy(cpl_plugin * plugin)
{
    cpl_recipe *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1 ;

    cpl_parameterlist_delete(recipe->parameters);
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Interpret the command line options and execute the data processing
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     oversample, swath_width, sum_only,
                            reduce_det, reduce_order, reduce_trace ;
    double                  smooth_slit ;
    cpl_frame           *   fr ;
    const char          *   science_file ;
    const char          *   trace_file ;

    hdrl_image          *   model_master[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slit_func_tab[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extract_tab[CR2RES_NB_DETECTORS] ;
    cpl_table           *   trace_table ;
    cpl_image           *   science_ima ;
    cpl_vector          **  spectrum ;
    cpl_vector          **  slit_func ;
    cpl_polynomial      **  traces ;
    hdrl_image          *   model_tmp ;
    cpl_vector          *   y_center ;
    int                     det_nr, extr_height, nb_traces, trace_id,
                            order, i ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.oversample");
    oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.swath_width");
    swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.smooth_slit");
    smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.sum_only");
    sum_only = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.height");
    extr_height = cpl_parameter_get_int(param);

    /* Check Parameters */
    /* TODO */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    fr = cpl_frameset_get_position(frameset, 0);
    science_file = cpl_frame_get_filename(fr) ;
    fr = cpl_frameset_get_position(frameset, 1);
    trace_file = cpl_frame_get_filename(fr) ;
    if (science_file == NULL || trace_file == NULL) {
        cpl_msg_error(__func__, "The utility needs a science file and a trace");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        model_master[det_nr-1] = NULL ;
        slit_func_tab[det_nr-1] = NULL ;
        extract_tab[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the trace table of this detector */
        cpl_msg_info(__func__, "Load the trace table") ;
        if ((trace_table = cr2res_io_load_TRACE_WAVE(trace_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,
                    "Failed to get trace table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        nb_traces = cpl_table_get_nrow(trace_table) ;

        /* Load the image in which the traces are to extract */
        if ((science_ima = cpl_image_load(science_file, CPL_TYPE_FLOAT,
                        0, det_nr)) == NULL) {
            cpl_table_delete(trace_table) ;
            cpl_msg_error(__func__, "Failed to load the image - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }

        /* Allocate Data containers */
        spectrum = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;
        slit_func = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;
        model_master[det_nr-1] = hdrl_image_create(science_ima, NULL) ;
        hdrl_image_mul_scalar(model_master[det_nr-1], (hdrl_value){0.0, 0.0}) ;

        /* Loop over the traces and extract them */
        for (i=0 ; i<nb_traces ; i++) {
            /* Initialise */
            slit_func[i] = NULL ;
            spectrum[i] = NULL ;
            model_tmp = NULL ;

            /* Get Order and trace id */
            order = cpl_table_get(trace_table, "Order", i, NULL) ;
            trace_id = cpl_table_get(trace_table, "TraceNb", i, NULL) ;

            /* Check if this order needs to be skipped */
            if (reduce_order > -1 && order != reduce_order) {
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Check if this trace needs to be skipped */
            if (reduce_trace > -1 && trace_id != reduce_trace) {
                cpl_msg_indent_less() ;
                continue ;
            }

            cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
            cpl_msg_indent_more() ;

            /* Get the 2 Edges for the current trace */
            if ((traces = cr2res_trace_open_get_polynomials(trace_table,
                            order, trace_id)) == NULL) {
                cpl_msg_warning(__func__, "Failed to get the traces");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Get the values between the 2 traces and the height */
            y_center = cr2res_trace_compute_middle(traces[0], traces[1],
                    cpl_image_get_size_x(science_ima)) ;
            if (extr_height == -1) {
                /* Only overwrite when input parameter not set*/
                extr_height = cr2res_trace_compute_height(traces[0], traces[1],
                        cpl_image_get_size_x(science_ima)) ;
            }
            cpl_polynomial_delete(traces[0]) ;
            cpl_polynomial_delete(traces[1]) ;
            cpl_free(traces) ;

            /* Call the SLIT DECOMPOSITION */
            if (cr2res_slitdec_vert(science_ima, y_center, extr_height,
                    swath_width, oversample, smooth_slit,
                    &(slit_func[i]), &(spectrum[i]), &model_tmp) != 0) {
                cpl_msg_error(__func__, "Cannot extract the trace") ;
                cpl_vector_delete(y_center) ;
                slit_func[i] = NULL ;
                spectrum[i] = NULL ;
                model_tmp = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            cpl_vector_delete(y_center) ;

            /* Update the model global image */
            if (model_tmp != NULL) {
                hdrl_image_add_image(model_master[det_nr-1], model_tmp) ;
                hdrl_image_delete(model_tmp) ;
            }
            cpl_msg_indent_less() ;
        }
        cpl_image_delete(science_ima) ;

        /* Create the slit_func_tab for the current detector */
        slit_func_tab[det_nr-1] = cr2res_slit_func_tab_create(
                slit_func, trace_table) ;

        /* Create the extracted_tab for the current detector */
        extract_tab[det_nr-1] = cr2res_extract_tab_create(
                spectrum, trace_table) ;
        cpl_table_delete(trace_table) ;

		/* Deallocate Vectors */
        for (i=0 ; i<nb_traces ; i++) {
            if (slit_func[i] != NULL) cpl_vector_delete(slit_func[i]) ;
            if (spectrum[i] != NULL) cpl_vector_delete(spectrum[i]) ;
        }
        cpl_free(spectrum) ;
        cpl_free(slit_func) ;
        cpl_msg_indent_less() ;
    }

    /* Save the Products */
    cr2res_io_save_SLIT_MODEL("cr2res_util_extract_model.fits", frameset,
            parlist, model_master, NULL, RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC("cr2res_util_extract_slit_func.fits", frameset,
            parlist, slit_func_tab, NULL, RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D("cr2res_util_extract_extract_1D.fits", frameset,
            parlist, extract_tab, NULL, RECIPE_STRING) ;

    /* Free and return */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        cpl_table_delete(slit_func_tab[det_nr-1]) ;
        cpl_table_delete(extract_tab[det_nr-1]) ;
        hdrl_image_delete(model_master[det_nr-1]) ;
    }
    return (int)cpl_error_get_code();
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the extract 1D table to be saved
  @param    spectrum   	The extracted spectra of the different orders
  @param    trace_table The trace table
  @return   the extract_1D table or NULL
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_extract_tab_create(
        cpl_vector      **  spectrum,
        cpl_table       *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pspec ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (spectrum == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check if all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL) {
            nrows = cpl_vector_get_size(spectrum[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL && cpl_vector_get_size(spectrum[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
		order = cpl_table_get(trace_table, "Order", i, NULL) ;
		trace_id = cpl_table_get(trace_table, "TraceNb", i, NULL) ;
        col_name = cpl_sprintf("%02d_%02d_SPEC", order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (spectrum[i] != NULL) {
            order = cpl_table_get(trace_table, "Order", i, NULL) ;
            trace_id = cpl_table_get(trace_table, "TraceNb", i, NULL) ;
            pspec = cpl_vector_get_data_const(spectrum[i]) ;
            col_name = cpl_sprintf("%02d_%02d_SPEC", order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pspec) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Create the slit functions table to be saved
  @param    slit_func   The slit functions of the different orders
  @param    trace_table The trace table
  @return   the slit_func table or NULL
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_slit_func_tab_create(
        cpl_vector      **  slit_func,
        cpl_table       *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pslit ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (slit_func == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check the all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (slit_func[i] != NULL) {
            nrows = cpl_vector_get_size(slit_func[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (slit_func[i] != NULL && cpl_vector_get_size(slit_func[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
		order = cpl_table_get(trace_table, "Order", i, NULL) ;
		trace_id = cpl_table_get(trace_table, "TraceNb", i, NULL) ;
        col_name = cpl_sprintf("%02d_%02d_SLIT_FUNC", order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (slit_func[i] != NULL) {
            order = cpl_table_get(trace_table, "Order", i, NULL) ;
            trace_id = cpl_table_get(trace_table, "TraceNb", i, NULL) ;
            pslit = cpl_vector_get_data_const(slit_func[i]) ;
            col_name = cpl_sprintf("%02d_%02d_SLIT_FUNC", order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pslit) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}
