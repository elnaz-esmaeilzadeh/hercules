/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "psolve.h"
#include "octor.h"
#include "util.h"
#include "stiffness.h"
#include "quake_util.h"
#include "cvm.h"
#include "drm_halfspace.h"
#include "topography.h"
#include "geometrics.h"


static int32_t            *myDRMFace1ElementsMapping;
static int32_t            *myDRMFace2ElementsMapping;
static int32_t            *myDRMFace3ElementsMapping;
static int32_t            *myDRMFace4ElementsMapping;
static int32_t            *myDRMBottomElementsMapping;

static int32_t            *myDRMBorder1ElementsMapping;
static int32_t            *myDRMBorder2ElementsMapping;
static int32_t            *myDRMBorder3ElementsMapping;
static int32_t            *myDRMBorder4ElementsMapping;

static int32_t            *myDRMBorder5ElementsMapping;
static int32_t            *myDRMBorder6ElementsMapping;
static int32_t            *myDRMBorder7ElementsMapping;
static int32_t            *myDRMBorder8ElementsMapping;

static int32_t            myTopoDRMFaceElementsCount = 0;
static int32_t            myTopoDRMBorderElementsCount = 0;


static double            myTs  = 0;
static double            myFc  = 0;
static int               myDir = 0;
static double            theDRMdepth;

static int32_t            myDRM_Face1Count  = 0;
static int32_t            myDRM_Face2Count  = 0;
static int32_t            myDRM_Face3Count  = 0;
static int32_t            myDRM_Face4Count  = 0;
static int32_t            myDRM_BottomCount = 0;
static int32_t            myDRM_Brd1 = 0;
static int32_t            myDRM_Brd2 = 0;
static int32_t            myDRM_Brd3 = 0;
static int32_t            myDRM_Brd4 = 0;
static int32_t            myDRM_Brd5 = 0;
static int32_t            myDRM_Brd6 = 0;
static int32_t            myDRM_Brd7 = 0;
static int32_t            myDRM_Brd8 = 0;

void drmHS_solver_init( mesh_t *myMesh, mysolver_t *mySolver) {


	int32_t halfFace_elem, theBorderElem;  /* half the number of elements in a box's face. Artificially defined by me */
	int32_t theFaceElem, theBaseElem;

	/* Data required by me before running the simmulation */
	double hmin     = 15.625 ; /* 15.625 for the SanchezSeama Ridge
	                              15.625 for the Gaussian Ridge */
	halfFace_elem   = 60;     /* this defines B. See figure in manuscript.
	                             halfFace_elem   = 15 for the SphericalRidge
	                             halfFace_elem   = 60 for the GaussianRidge
	                             halfFace_elem   = 60 for the SanSesma Ridge */

	theBorderElem   = 5;     /* this defines the depth, See figure in manuscript  */
	                         /* 5 for the gaussian ridge, 3 for the SanSesma ridge */
	double theXc           = 1000;   /* domain center coordinates */
	double theYc           = 1000;
	double Ts              = 0.20;  /* 0.18 for the Gaussian Ridge. 0.2 for the SanSesma Ridge  */
	double fc              = 10;  /* 10.26 for the Gaussian Ridge. 10 for the SanSesma Ridge  */
	int    wave_dir        = 2;    /*  0:X, 1:Y, 2:Z */

	/* --------  */
	double DRM_D = theBorderElem * hmin;
	double DRM_B = halfFace_elem * hmin;

	double thebase_zcoord = get_thebase_topo();

	theFaceElem = 2 * ( halfFace_elem + 0 ) * theBorderElem;
	theBaseElem = 4 * halfFace_elem * halfFace_elem;
	theDRMdepth = DRM_D;

	myTopoDRMBorderElementsCount = theBorderElem;
	myTopoDRMFaceElementsCount   = theFaceElem;
	myTs                         = Ts;
	myFc                         = fc;
	myDir                        = wave_dir;


	/*  mapping of face1 elements */
	int32_t eindex;
	int32_t countf1 = 0, countf2 = 0, countf3 = 0, countf4 = 0, countbott=0;
	int32_t countb1 = 0, countb2 = 0, countb3 = 0, countb4 = 0;
	int32_t countb5 = 0, countb6 = 0, countb7 = 0, countb8 = 0;

	XMALLOC_VAR_N(myDRMFace1ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace2ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace3ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMFace4ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBottomElementsMapping , int32_t, theBaseElem);

	/* border elements*/
	XMALLOC_VAR_N(myDRMBorder1ElementsMapping, int32_t, theBorderElem + 1);
	XMALLOC_VAR_N(myDRMBorder2ElementsMapping, int32_t, theBorderElem + 1);
	XMALLOC_VAR_N(myDRMBorder3ElementsMapping, int32_t, theBorderElem + 1);
	XMALLOC_VAR_N(myDRMBorder4ElementsMapping, int32_t, theBorderElem + 1);

	XMALLOC_VAR_N(myDRMBorder5ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder6ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder7ElementsMapping, int32_t, theFaceElem);
	XMALLOC_VAR_N(myDRMBorder8ElementsMapping, int32_t, theFaceElem);

	for (eindex = 0; eindex < myMesh->lenum; eindex++) {

		elem_t     *elemp;
		node_t     *node0dat;
		edata_t    *edata;
		double      xo, yo, zo;
		int32_t	    node0;

		elemp    = &myMesh->elemTable[eindex]; //Takes the information of the "eindex" element
		node0    = elemp->lnid[0];             //Takes the ID for the zero node in element eindex
		node0dat = &myMesh->nodeTable[node0];
		edata    = (edata_t *)elemp->data;

		/* get coordinates of element zero node */
		xo = (node0dat->x)*(myMesh->ticksize);
		yo = (node0dat->y)*(myMesh->ticksize);
		zo = (node0dat->z)*(myMesh->ticksize);


		if ( ( ( yo - theYc ) == DRM_B ) &&                            /*face 1*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace1ElementsMapping[countf1] = eindex;
			countf1++;
		} else 	if ( ( ( theYc - ( yo + hmin ) ) == DRM_B ) &&        /* face 2*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace2ElementsMapping[countf2] = eindex;
			countf2++;

		} else 	if ( ( ( theXc - ( xo + hmin ) ) == DRM_B ) &&      /*face 3*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace3ElementsMapping[countf3] = eindex;
			countf3++;

		} else 	if ( ( ( xo - theXc ) == DRM_B ) &&             /* face 4*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo <  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMFace4ElementsMapping[countf4] = eindex;
			countf4++;

		} else 	if ( ( yo >= ( theYc - DRM_B ) ) &&         /*bottom*/
				( yo <  ( theYc + DRM_B ) ) &&
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBottomElementsMapping[countbott] = eindex;
			countbott++;

		} else 	if ( ( ( yo - theYc ) == DRM_B ) &&      /*border 1*/
				( ( xo + hmin ) == ( theXc - DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder1ElementsMapping[countb1] = eindex;
			countb1++;

		} else if ( ( ( yo - theYc ) == DRM_B  ) &&      /*border 2*/
				( xo  == ( theXc + DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder2ElementsMapping[countb2] = eindex;
			countb2++;

		} else if ( ( ( theYc - ( yo + hmin ) ) == DRM_B ) &&  /* border 3*/
				( ( xo + hmin)  == ( theXc - DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder3ElementsMapping[countb3] = eindex;
			countb3++;

		} else if ( ( ( theYc - ( yo + hmin ) ) == DRM_B ) &&  /* border 4*/
				(  xo  == ( theXc + DRM_B ) ) &&
				( zo <=  DRM_D + thebase_zcoord ) &&
				( zo >=  thebase_zcoord ) ) {

			myDRMBorder4ElementsMapping[countb4] = eindex;
			countb4++;

		} else 	if ( ( ( yo - theYc ) == DRM_B ) &&            /* border 5*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder5ElementsMapping[countb5] = eindex;
			countb5++;

		} else if ( ( ( theYc - ( yo + hmin ) ) == DRM_B ) &&          /* border 6*/
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder6ElementsMapping[countb6] = eindex;
			countb6++;

		} else if ( ( ( theXc - ( xo + hmin ) ) == DRM_B ) &&      /* border 7*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder7ElementsMapping[countb7] = eindex;
			countb7++;

		} else if ( ( ( xo - theXc ) == DRM_B ) &&             /* border 8*/
				( yo >= ( theYc - DRM_B ) ) &&
				( yo <  ( theYc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBorder8ElementsMapping[countb8] = eindex;
			countb8++;
		}
	}


	myDRM_Face1Count  = countf1;
	myDRM_Face2Count  = countf2;
	myDRM_Face3Count  = countf3;
	myDRM_Face4Count  = countf4;
	myDRM_BottomCount = countbott;
	myDRM_Brd1        = countb1;
	myDRM_Brd2        = countb2;
	myDRM_Brd3        = countb3;
	myDRM_Brd4        = countb4;
	myDRM_Brd5        = countb5;
	myDRM_Brd6        = countb6;
	myDRM_Brd7        = countb7;
	myDRM_Brd8        = countb8;

}

