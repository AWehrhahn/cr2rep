/*
 * This file is part of the CR2RE Pipeline
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
#include "cr2res_slitdec.h"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_extract_create(cpl_plugin *);
static int cr2res_util_extract_exec(cpl_plugin *);
static int cr2res_util_extract_destroy(cpl_plugin *);
static int cr2res_util_extract(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_extract_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"science.fits " CR2RES_OBS_1D "\n"
"trace.fits " CR2RES_TRACE_OPEN_PROCATG "\n"
"trace.fits " CR2RES_TRACE_DECKER_PROCATG "\n"
"trace-decker.fits " CR2RES_FLAT_DECKER_RAW "\n"
"master_bpm.fits " CR2RES_MASTER_BPM_PROCATG "\n"
" The recipe produces the following products:\n"
"\n";

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
                    "cr2res_util_extract",
                    "extract utility recipe",
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
    p = cpl_parameter_new_value("cr2res.cr2res_extract.oversample",
            CPL_TYPE_INT,
            "factor by which to oversample the extraction",
            "cr2res.cr2res_extract", 10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_extract.swath_width",
            CPL_TYPE_INT,
            "The swath width (number of columns) over which the slit function is assumed constant",
            "cr2res.cr2res_extract", 256);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "swathwidth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_extract.sum_only",
            CPL_TYPE_BOOL,
            "If True, sum along detector column only, instead of slit decomposition",
            "cr2res.cr2res_extract", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth");
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
    cpl_propertylist    *   plist;
    cpl_frame           *   rawframe ;
    cpl_propertylist    *   applist;
    cpl_image           *   in ;
    cpl_image           *   trace ;
    const cpl_parameter * param;
    cpl_frameset        *   openslit_frames;
    cpl_frameset        *   decker_frames;
    int                     npolys;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.min_cluster");
    int mincluster = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_extract.poly_order");
    int polyorder = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_extract.smooth");
    int smoothfactor = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_extract.cpl-lab");
    int cpl_lab = cpl_parameter_get_bool(param);

    /* Check Parameters */
    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Data */
    openslit_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_OPEN_RAW);
    int nb_open = cpl_frameset_get_size(openslit_frames);
    decker_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_DECKER_RAW);
    int nb_decker = cpl_frameset_get_size(decker_frames);
    cpl_msg_info(__func__, "Got %d & %d open slit and decker files, resp.", nb_open, nb_decker);

    int i;
    for (i=0; i<nb_open; i++){
        rawframe = cpl_frameset_get_position(frameset, i);
        in = cpl_image_load(cpl_frame_get_filename(rawframe), CPL_TYPE_DOUBLE,0,0);
        if (in == NULL) {
            cpl_msg_error(__func__, "Cannot load the input image") ;
            cpl_propertylist_delete(plist) ;
            return -1 ;
        }
        plist = cpl_propertylist_load(cpl_frame_get_filename(rawframe), 0);
        model = cr2res_slitdec_vert(in, CR2RES_DECKER_NONE, &npolys); //TODO: fix call
        cpl_msg_debug(__func__,"%d npolys found.",npolys);
    }
    cpl_frameset_delete(openslit_frames);
    for (i=0; i<nb_decker; i++){
        rawframe = cpl_frameset_get_position(frameset, i);
    }
    cpl_frameset_delete(decker_frames);


    trace = cpl_image_duplicate(in);
    if (trace == NULL) {
        cpl_propertylist_delete(plist) ;
        cpl_image_delete(in) ;
        return -1 ;
    }


    cpl_image_delete(in) ;

    /* Add the product category  */
    applist = cpl_propertylist_duplicate(plist);
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG,
            CR2RES_TRACE_OPEN_PROCATG);

    /* Save Product */
    cpl_dfs_save_image(frameset, plist, parlist, frameset, NULL,
            trace, CPL_BPP_IEEE_FLOAT, "cr2res_util_extract", applist,
            NULL, PACKAGE "/" PACKAGE_VERSION, "cr2res_util_extract.fits") ;

    /* Free and return */
    cpl_propertylist_delete(plist) ;
    cpl_propertylist_delete(applist) ;
    cpl_image_delete(trace) ;
    return (int)cpl_error_get_code();
}
