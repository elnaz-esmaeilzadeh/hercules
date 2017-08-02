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

static planewavetype_t	thePlaneWaveType;
static int32_t	        theDRMBox_halfwidthElements = 0;
static int32_t	        theDRMBox_DepthElements = 0;
static double 	        thedrmbox_esize         = 0.0;

static double 	        theTs = 0.0;
static double 	        thefc = 0.0;
static double           theUo = 0.0;
static double 	        theplanewave_strike = 0.0;
static double 	        theXc  = 0.0;
static double 	        theYc  = 0.0;


static int32_t            *myDRMBottomElementsMapping;
static double              theDRMdepth;
static int32_t            myDRM_BottomCount = 0;


void drm_planewaves_init ( int32_t myID, const char *parametersin ) {

    int     int_message[3];
    double  double_message[7];

    /* Capturing data from file --- only done by PE0 */
    if (myID == 0) {
        if ( drm_planewaves_initparameters( parametersin ) != 0 ) {
            fprintf(stderr,"Thread %d: drm_planewaves_init: "
                    "incidentPlaneWaves_initparameters error\n",myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }

    }

    /* Broadcasting data */

    int_message   [0]    = (int)thePlaneWaveType;
    int_message   [1]    = theDRMBox_halfwidthElements;
    int_message   [2]    = theDRMBox_DepthElements;

    double_message[0] = theTs;
    double_message[1] = thefc;
    double_message[2] = theUo;
    double_message[3] = theplanewave_strike;
    double_message[4] = theXc;
    double_message[5] = theYc;
    double_message[6] = thedrmbox_esize;

    MPI_Bcast(double_message, 7, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(int_message,    3, MPI_INT,    0, comm_solver);

    thePlaneWaveType             = int_message[0];
    theDRMBox_halfwidthElements  = int_message[1];
    theDRMBox_DepthElements      = int_message[2];

    theTs               = double_message[0];
    thefc               = double_message[1];
    theUo               = double_message[2];
    theplanewave_strike = double_message[3];
    theXc               = double_message[4];
    theYc               = double_message[5];
    thedrmbox_esize     = double_message[6];

    return;

}



int32_t
drm_planewaves_initparameters ( const char *parametersin ) {
    FILE                *fp;

    double              drmbox_halfwidth_elements, drmbox_depth_elements, Ts, fc, Uo, planewave_strike, L_ew, L_ns, drmbox_esize;
    char                type_of_wave[64];

    planewavetype_t     planewave;


    /* Opens parametersin file */

   if ( ( fp = fopen(parametersin, "r" ) ) == NULL ) {
        fprintf( stderr,
                 "Error opening %s\n at drm_planewaves_initparameters",
                 parametersin );
        return -1;
    }


     /* Parses parametersin to capture drm_planewaves single-value parameters */
    if ( ( parsetext(fp, "type_of_wave",                     's', &type_of_wave                ) != 0) ||
         ( parsetext(fp, "DRMBox_Noelements_in_halfwidth",   'd', &drmbox_halfwidth_elements   ) != 0) ||
         ( parsetext(fp, "DRMBox_Noelements_in_depth",       'd', &drmbox_depth_elements       ) != 0) ||
         ( parsetext(fp, "DRMBox_element_size_m",            'd', &drmbox_esize                ) != 0) ||
         ( parsetext(fp, "Ts",                               'd', &Ts                          ) != 0) ||
         ( parsetext(fp, "region_length_east_m",             'd', &L_ew                        ) != 0) ||
         ( parsetext(fp, "region_length_north_m",            'd', &L_ns                        ) != 0) ||
         ( parsetext(fp, "fc",                               'd', &fc                          ) != 0) ||
         ( parsetext(fp, "Uo",                               'd', &Uo                          ) != 0) ||
         ( parsetext(fp, "planewave_strike",                 'd', &planewave_strike            ) != 0) )
    {
        fprintf( stderr,
                 "Error parsing planewaves parameters from %s\n",
                 parametersin );
        return -1;
    }

    if ( strcasecmp(type_of_wave, "SV") == 0 ) {
    	planewave = SV;
    } else if ( strcasecmp(type_of_wave, "P") == 0 ) {
    	planewave = P;
    } else {
        fprintf(stderr,
                "Illegal type_of_wave for incident plane wave analysis"
                "(SV, P): %s\n", type_of_wave);
        return -1;
    }

    /*  Initialize the static global variables */
	thePlaneWaveType                 = planewave;
	theDRMBox_halfwidthElements      = drmbox_halfwidth_elements;
	theDRMBox_DepthElements          = drmbox_depth_elements;
	theTs                            = Ts;
	thefc                            = fc;
    theUo                            = Uo;
	theplanewave_strike              = planewave_strike * PI / 180.00;
	theXc                            = L_ew / 2.0;
	theYc                            = L_ns / 2.0;
	thedrmbox_esize                  = drmbox_esize;

    fclose(fp);

    return 0;
}


void PlaneWaves_solver_init( int32_t myID, mesh_t *myMesh, mysolver_t *mySolver) {

	int32_t theFaceElem, theBaseElem;

	double DRM_D = theDRMBox_DepthElements * thedrmbox_esize;
	double DRM_B = theDRMBox_halfwidthElements * thedrmbox_esize;
	double thebase_zcoord = get_thebase_topo();

	theFaceElem = 2 * ( theDRMBox_halfwidthElements + 0 ) * theDRMBox_DepthElements;
	theBaseElem = 4 * theDRMBox_halfwidthElements * theDRMBox_halfwidthElements;
	theDRMdepth = DRM_D;

	/*  mapping of face1 elements */
	int32_t eindex;
	int32_t countbott=0;

	XMALLOC_VAR_N(myDRMBottomElementsMapping , int32_t, theBaseElem);

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

		if ( ( yo >= ( theYc - DRM_B ) ) &&         /*bottom*/
				( yo <  ( theYc + DRM_B ) ) &&
				( xo >= ( theXc - DRM_B ) ) &&
				( xo <  ( theXc + DRM_B ) ) &&
				( zo ==  DRM_D + thebase_zcoord )          ) {

			myDRMBottomElementsMapping[countbott] = eindex;
			countbott++;

		}

	}

	myDRM_BottomCount = countbott;
}


void compute_addforce_PlaneWaves ( mesh_t     *myMesh,
                                mysolver_t *mySolver,
                                double      theDeltaT,
                                int         step,
                                fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8])
{

    int32_t   eindex;
    int32_t   face_eindex;

    int  f_nodes_bottom[4] = { 0, 1, 2, 3 };
    int  e_nodes_bottom[4] = { 4, 5, 6, 7 };

    double tt = theDeltaT * step;
    int  *f_nodes, *e_nodes;

    /* Loop over bottom elements */
    f_nodes = &f_nodes_bottom[0];
    e_nodes = &e_nodes_bottom[0];
    for ( face_eindex = 0; face_eindex < myDRM_BottomCount ; face_eindex++) {
    	eindex = myDRMBottomElementsMapping[face_eindex];
    	DRM_ForcesinElement ( myMesh, mySolver, theK1, theK2, f_nodes, e_nodes, eindex, tt, 4, 4 );
    } /* all elements in bottom*/


    return;
}


void DRM_ForcesinElement ( mesh_t     *myMesh,
		mysolver_t *mySolver,
		fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8],
		int *f_nodes, int *e_nodes, int32_t   eindex, double tt, int Nnodes_e, int Nnodes_f )
{

    int       i, j;
    int  CoordArr[8]      = { 0, 0, 0, 0, 1, 1, 1, 1 };

    double thebase_zcoord = get_thebase_topo();

	fvector_t localForce[8];

	elem_t        *elemp;
	edata_t       *edata;
	node_t        *node0dat;
	double        xo, yo, zo;
	int32_t	      node0;
	e_t*          ep;

	/* Capture the table of elements from the mesh and the size
	 * This is what gives me the connectivity to nodes */
	elemp        = &myMesh->elemTable[eindex];
	edata        = (edata_t *)elemp->data;
	node0        = elemp->lnid[0];
	node0dat     = &myMesh->nodeTable[node0];
	ep           = &mySolver->eTable[eindex];

	/* get coordinates of element zero node */
	xo = (node0dat->x)*(myMesh->ticksize);
	yo = (node0dat->y)*(myMesh->ticksize);
	zo = (node0dat->z)*(myMesh->ticksize) - thebase_zcoord;


	/* get material properties  */
	double  h, Vs;
	h    = (double)edata->edgesize;

	if ( thePlaneWaveType == SV  )
		Vs = edata->Vs;
	else
		Vs = edata->Vp;


	/* Force contribution from external nodes */
	/* -------------------------------
	 * Ku DONE IN THE CONVENTIONAL WAY
	 * ------------------------------- */
	memset( localForce, 0, 8 * sizeof(fvector_t) );

	fvector_t myDisp;
	/* forces over f nodes */
	for (i = 0; i < Nnodes_f; i++) {

		int  nodef = *(f_nodes + i);
		fvector_t* toForce = &localForce[ nodef ];

		/* incoming displacements over e nodes */
		for (j = 0; j < Nnodes_e; j++) {

			int  nodee = *(e_nodes + j);

			double z_ne = zo + h * CoordArr[ nodee ];   /* get zcoord */
			getRicker ( &myDisp, z_ne, tt, Vs ); /* get Displ */

			MultAddMatVec( &theK1[ nodef ][ nodee ], &myDisp, -ep->c1, toForce );
			MultAddMatVec( &theK2[ nodef ][ nodee ], &myDisp, -ep->c2, toForce );
		}
	}

	/* forces over e nodes */
	for (i = 0; i < Nnodes_e; i++) {

		int  nodee = *(e_nodes + i);
		fvector_t* toForce = &localForce[ nodee ];

		/* incoming displacements over f nodes */
		for (j = 0; j < Nnodes_f; j++) {

			int  nodef = *(f_nodes + j);

			double z_nf = zo + h * CoordArr[ nodef ];   /* get zcoord */
			getRicker ( &myDisp, z_nf, tt, Vs ); /* get Displ */

			MultAddMatVec( &theK1[ nodee ][ nodef ], &myDisp, ep->c1, toForce );
			MultAddMatVec( &theK2[ nodee ][ nodef ], &myDisp, ep->c2, toForce );
		}
	}
	/* end Ku */

	/* Loop over the 8 element nodes:
	 * Add the contribution calculated above to the node
	 * forces carried from the source and stiffness.
	 */
	for (i = 0; i < 8; i++) {

		int32_t    lnid;
		fvector_t *nodalForce;

		lnid = elemp->lnid[i];

		nodalForce = mySolver->force + lnid;

		nodalForce->f[0] += localForce[i].f[0] ;
		nodalForce->f[1] += localForce[i].f[1] ;
		nodalForce->f[2] += localForce[i].f[2] ;

	} /* element nodes */

}

void getRicker ( fvector_t *myDisp, double zp, double t, double Vs ) {

	double Rz = Ricker_displ ( zp, theTs, t, thefc, Vs  );

	if ( thePlaneWaveType == SV ) {
		myDisp->f[0] = Rz * theUo * cos (theplanewave_strike);
		myDisp->f[1] = Rz * theUo * sin (theplanewave_strike);
		myDisp->f[2] = 0.0;
	} else {
		myDisp->f[0] = 0.0;
		myDisp->f[1] = 0.0;
		myDisp->f[2] = Rz * theUo;
	}

}

double Ricker_displ ( double zp, double Ts, double t, double fc, double Vs  ) {

	double alfa1 = ( PI * fc ) * ( PI * fc ) * ( t - zp / Vs - Ts) * ( t - zp / Vs - Ts);
	double alfa2 = ( PI * fc ) * ( PI * fc ) * ( t + zp / Vs - Ts) * ( t + zp / Vs - Ts);

	double uo1 = ( 2.0 * alfa1 - 1.0 ) * exp(-alfa1);
	double uo2 = ( 2.0 * alfa2 - 1.0 ) * exp(-alfa2);

	return (uo1+uo2);
}

