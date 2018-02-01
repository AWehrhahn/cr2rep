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

#ifndef CR2RES_IO_H
#define CR2RES_IO_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

int cr2res_io_get_ext_idx(
        const char  *   filename,
        int             detector) ;

cpl_bivector * cr2res_io_load_EMISSION_LINES(
        const char  *   filename) ;

cpl_image * cr2res_io_load_MASTER_DARK(
        const char  *   filename,
        int             detector,
        int             data);

cpl_image * cr2res_io_load_MASTER_BPM(
        const char  *   filename,
        int             detector);

cpl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_MASTER_FLAT(
        const char  *   filename,
        int             detector,
        int             data);

cpl_table * cr2res_io_load_TRACE_WAVE(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_BLAZE(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_BLAZE_IMAGE(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_SLIT_MODEL(
        const char  *   filename,
        int             detector,
        int             data);

cpl_table * cr2res_io_load_SLIT_FUNC(
        const char  *   filename,
        int             detector) ;

cpl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_SLITPOS_MAP(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_TILT_MAP(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_TILT_POLY(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_SPLICED_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_2D(
        const char  *   filename);

cpl_table * cr2res_io_load_EXTRACT_POL(
        const char  *   filename,
        int             detector);

int cr2res_io_save_EMISSION_LINES(
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set,
        const char              *   recipe) ;

int cr2res_io_save_MASTER_DARK(
        cpl_frameset            *   allframes,
        const char              *   filename,
        cpl_frameset            *   used_frames,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe) ;

int cr2res_io_save_MASTER_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_bpms,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_imagelist          **  detlin_coeffs,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_MASTER_FLAT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_flats,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_TRACE_WAVE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_BLAZE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_BLAZE_IMAGE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  blaze,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_SLIT_MODEL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_SLIT_FUNC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_SLIT_ILLUM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe) ;

int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_SLITPOS_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_TILT_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_TILT_POLY(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_EXTRACT_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe);

int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   table,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe);

int cr2res_io_save_EXTRACT_POL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe) ;

#endif

