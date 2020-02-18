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

#include "cr2res_utils.h"
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"
#include "cr2res_qc.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_obs_staring"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_obs_staring_find_RAW(
        const cpl_frameset  *   in) ;
static int cr2res_obs_staring_check_inputs_validity(
        const cpl_frameset  *   rawframes) ;
static int cr2res_obs_staring_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     reduce_det,
        cpl_table           **  extract,
        cpl_table           **  slitfunc,
        hdrl_image          **  model,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_obs_staring_create(cpl_plugin *);
static int cr2res_obs_staring_exec(cpl_plugin *);
static int cr2res_obs_staring_destroy(cpl_plugin *);
static int cr2res_obs_staring(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_staring_description[] = "\
Staring Observation                                                     \n\
  This recipe handles staring observations.                             \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_STARING_OTHER_RAW" [1 to n]                   \n\
          or " CR2RES_OBS_STARING_JITTER_RAW" [1 to n]                  \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1]                       \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    detlin.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG " [0 to 1]           \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_MERGE_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    master_flat.fits " CR2RES_CAL_FLAT_MASTER_PROCATG " [0 to 1]        \n\
                                                                        \n\
  Outputs   TODO                                                        \n\
                                                                        \n\
  Algorithm      TODO                                                   \n\
                                                                        \n\
  Library functions uѕed   TODO                                         \n\
";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

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
                    RECIPE_STRING,
                    "Staring Observation recipe",
                    cr2res_obs_staring_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_staring_create,
                    cr2res_obs_staring_exec,
                    cr2res_obs_staring_destroy)) {    
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
static int cr2res_obs_staring_create(cpl_plugin * plugin)
{
    cpl_recipe          *   recipe ;
    cpl_parameter       *   p ;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_staring", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_staring", 90);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_obs_staring", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_smooth",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_obs_staring", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_staring", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
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
static int cr2res_obs_staring_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_staring(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_staring(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     extract_oversample, extract_swath_width,
                            extract_height, reduce_det, ndit, nexp,
                            disp_order, disp_trace ;
    double                  extract_smooth, ra, dec, dit ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   photo_flux_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    cpl_table           *   extract[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfunc[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   model[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   plist ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr, type; 


    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_smooth");
    extract_smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.display_order");
    disp_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.display_trace");
    disp_trace = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */
    trace_wave_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    if (trace_wave_frame == NULL) {
        cpl_msg_error(__func__, "Could not find TRACE_WAVE frame") ;
        return -1 ;
    }
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DARK_MASTER_PROCATG) ; 
    photo_flux_frame = cpl_frameset_find_const(frameset,
            CR2RES_PHOTO_FLUX_PROCATG) ; 
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_MASTER_PROCATG) ; 
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Get the RAW Frames */
    rawframes = cr2res_obs_staring_find_RAW(frameset) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
      
    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        extract[det_nr-1] = NULL ;
        slitfunc[det_nr-1] = NULL ;
        model[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_staring_reduce(rawframes, 
                    trace_wave_frame, detlin_frame, master_dark_frame, 
                    master_flat_frame, bpm_frame, 0, extract_oversample, 
                    extract_swath_width, extract_height, extract_smooth, det_nr,
                    &(extract[det_nr-1]),
                    &(slitfunc[det_nr-1]),
                    &(model[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    out_file = cpl_sprintf("%s_slitfunc.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC(out_file, frameset, rawframes, parlist,
            slitfunc, NULL, ext_plist, CR2RES_OBS_STARING_SLITFUNC_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_model.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_MODEL(out_file, frameset, rawframes, parlist,
            model, NULL, ext_plist, CR2RES_OBS_STARING_SLITMODEL_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_extracted.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, extract,
            NULL, ext_plist, CR2RES_OBS_STARING_EXTRACT_PROCATG, RECIPE_STRING);
    cpl_free(out_file);

    /* Free */
    cpl_frameset_delete(rawframes) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (extract[det_nr-1] != NULL) 
            cpl_table_delete(extract[det_nr-1]) ;
        if (slitfunc[det_nr-1] != NULL) 
            cpl_table_delete(slitfunc[det_nr-1]) ;
        if (model[det_nr-1] != NULL)
            hdrl_image_delete(model[det_nr-1]) ;
        if (ext_plist[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plist[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the science recipe on a specific detector
  @param rawframes              Raw science frames
  @param trace_wave_frame       Trace Wave file
  @param detlin_frame           Associated detlin coefficients
  @param master_dark_frame      Associated master dark
  @param master_flat_frame      Associated master flat
  @param bpm_frame              Associated BPM
  @param calib_cosmics_corr     Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth         Extraction related
  @param reduce_det             The detector to compute
  @param extract                [out] extracted spectrum 
  @param slitfunc               [out] slit function
  @param model                  [out] slit model
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     reduce_det,
        cpl_table           **  extract,
        cpl_table           **  slitfunc,
        hdrl_image          **  model,
        cpl_propertylist    **  ext_plist)
{
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    cpl_vector          *   dits ;
    cpl_table           *   trace_wave ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    cpl_propertylist    *   plist ;
    cpl_size                nframes, i ;
    hdrl_image          *   model_master ;
    cpl_table           *   slit_func ;
    cpl_table           *   extracted ;
    int                     det_nr ;

    /* Check Inputs */
    if (extract == NULL || ext_plist == NULL || rawframes == NULL
            || trace_wave_frame == NULL) return -1 ;

    /* Check raw frames consistency */
    if (cr2res_obs_staring_check_inputs_validity(rawframes) != 1) {
        cpl_msg_error(__func__, "Invalid Inputs") ;
        return -1 ;
    }

    /* Initialise */
    nframes = cpl_frameset_get_size(rawframes) ;

    /* Load the DITs if necessary */
    if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /* Load image list */
    if ((in = cr2res_io_load_image_list_from_set(rawframes, 
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0, 0, 
                    master_flat_frame, master_dark_frame, bpm_frame, 
                    detlin_frame, dits)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;

    /* Collapse the image list */
    cpl_msg_info(__func__, "Collapse") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(in_calib, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        hdrl_imagelist_delete(in_calib) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_image_delete(contrib) ;
    hdrl_imagelist_delete(in_calib) ;
    cpl_msg_indent_less() ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        hdrl_image_delete(collapsed) ;
        return -1 ;
    }

    /* Execute the extraction */
    cpl_msg_info(__func__, "Spectra Extraction") ;
    if (cr2res_extract_traces(collapsed, trace_wave, -1, -1,
                CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth,
                &extracted, &slit_func, &model_master) == -1) {
        cpl_msg_error(__func__, "Failed to extract");
        hdrl_image_delete(collapsed) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_table_delete(trace_wave) ;
	hdrl_image_delete(collapsed) ;

    /* QC parameters */
    plist = cpl_propertylist_new() ;

    /* Compute the QC parameters */

    /* Store the QC parameters in the plist */

    /* Return */
    *extract = extracted ;
    *slitfunc = slit_func ;
    *model = model_master ;
    *ext_plist = plist ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Run basic checks for the rawframes consistency
  @param    rawframes   The input rawframes
  @return   1 if valid, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_check_inputs_validity(
        const cpl_frameset  *   rawframes)
{
    cpl_propertylist        *   plist ;
    cr2res_nodding_pos      *   nod_positions ;
    cpl_size                    nframes, i ;
    double                      nodthrow, nodthrow_cur ;
    int                         nb_a, nb_b ;

    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Need even number of frames */
    nframes = cpl_frameset_get_size(rawframes) ;
    if (nframes % 2) {
        cpl_msg_error(__func__, "Require an even number of raw frames") ;
        return 0 ;
    }

    /* Need same number of A and B positions */
    nb_a = nb_b = 0 ;
    nod_positions = cr2res_nodding_read_positions(rawframes) ;
    if (nod_positions == NULL) return -1 ;
    for (i=0 ; i<nframes ; i++) {
        if (nod_positions[i] == CR2RES_NODDING_A) nb_a++ ;
        if (nod_positions[i] == CR2RES_NODDING_B) nb_b++ ;
    }
    cpl_free(nod_positions) ;    

    if (nb_a == 0 || nb_b == 0 || nb_a != nb_b) {
        cpl_msg_error(__func__, "Require as many A and B positions") ;
        return 0 ;
    }

    /* Need the same nod throw in all frames */
    if ((plist = cpl_propertylist_load(cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes, 0)),
                    0)) == NULL) {
        return -1;
    } 
    nodthrow = cr2res_pfits_get_nodthrow(plist);
    cpl_propertylist_delete(plist) ;
    for (i=1 ; i<nframes ; i++) {
        if ((plist = cpl_propertylist_load(cpl_frame_get_filename(
                            cpl_frameset_get_position_const(rawframes, i)),
                        0)) == NULL) {
            return -1;
        } 
        nodthrow_cur = cr2res_pfits_get_nodthrow(plist);
        cpl_propertylist_delete(plist) ;

        if (fabs(nodthrow_cur-nodthrow) > 1e-3) {
            cpl_msg_error(__func__, 
                    "Require constant NOD THROW in the raw frames") ;
            return 0 ;
        }
    }
    return 1 ;
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_OBS_STARING_OTHER_RAW
                        CR2RES_OBS_STARING_JITTER_RAW
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_obs_staring_find_RAW(
        const cpl_frameset  *   in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_OBS_STARING_OTHER_RAW) ;
    if (out == NULL) {
        out = cr2res_extract_frameset(in, CR2RES_OBS_STARING_JITTER_RAW) ;
    }
    return out ;
}
