#include <iostream>
#include "energy.h"
#include "graph.h"
#include "GCoptimization.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits>
#define MAX_INTT 1000000000

// <!--Msun
#ifdef MEX_COMPILE
#include "mex.h"
#endif

#include "QPBO.h"
// Msun -->

/**************************************************************************************/
// <!--Msun
void GCoptimization::UseQPBO(){
	m_solver = QPBOsolver;
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("USEQPBO\n");
#endif
#endif
}

void GCoptimization::SetSegWeight(EnergyTermType seg_weight){
	m_seg_weight = seg_weight;
}

void GCoptimization::SetDetClassWeights( double* det_class_weights){
	for (int i=0; i<m_num_labels; i++)
		m_det_class_weights[i] = (EnergyTermType)det_class_weights[i];
}
// Msun -->

/**************************************************************************************/

void GCoptimization::initialize_memory()
{
    int i = 0;
    
	m_lookupPixVar = (PixelType *) new PixelType[m_num_pixels];
	m_labelTable   = (LabelType *) new LabelType[m_num_labels];

	terminateOnError( !m_lookupPixVar || !m_labelTable,"Not enough memory");

	for ( i = 0; i < m_num_labels; i++ )
		m_labelTable[i] = i;

	for ( i = 0; i < m_num_pixels; i++ )
		m_lookupPixVar[i] = -1;

	// <!--Msun	
	m_lookupValidPix = (PixelType *) new PixelType[m_num_pixels];
	for ( i = 0; i < m_num_pixels; i++ )
		m_lookupValidPix[i] = -1;
	// Msun -->
}


/**************************************************************************************/

void GCoptimization::commonGridInitialization(PixelType width, PixelType height, int nLabels)
{
	terminateOnError( (width < 0) || (height <0) || (nLabels <0 ),"Illegal negative parameters");

	m_width              = width;
	m_height             = height;
	m_num_pixels         = width*height;
	m_num_labels         = nLabels;
	m_grid_graph         = 1;
		
	
	//srand(time(NULL));
	srand(0); //Msun: 3/25 fix randomness
}
/**************************************************************************************/

void GCoptimization::commonNonGridInitialization(PixelType nupixels, int num_labels)
{
	terminateOnError( (nupixels <0) || (num_labels <0 ),"Illegal negative parameters");

	m_num_labels		 = num_labels;
	m_num_pixels		 = nupixels;
	m_grid_graph         = 0;

#ifdef LOG_ON
	std::cout<<"m_num_pixels="<<m_num_pixels<<std::endl;
	std::cout<<"m_neighbors step to nupixels="<<nupixels<<std::endl;
#endif
	m_neighbors = (LinkedBlockList *) new LinkedBlockList[nupixels];

	terminateOnError(!m_neighbors,"Not enough memory");

}

/**************************************************************************************/

void GCoptimization::commonInitialization(int dataSetup, int smoothSetup)
{
	int i;

	m_random_label_order = 1;
	m_random_pair_thing_order = false; // should be always false
	m_dataInput          = dataSetup;
	m_smoothInput        = smoothSetup;


	//Min hack
	// Msun: 3/5 added loss
	m_losscost    = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_pixels];
	for ( i = 0; i < m_num_labels*m_num_pixels; i++ ) 
		m_losscost[i] = (EnergyTermType) 0;

	if (m_dataInput == SET_INDIVIDUALLY )
	{
		m_datacost    = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_pixels];
		terminateOnError(!m_datacost,"Not enough memory");
		for ( i = 0; i < m_num_labels*m_num_pixels; i++ ) 
			m_datacost[i] = (EnergyTermType) 0;
		m_dataType = ARRAY;
	}
	else m_dataType = NONE;


	if ( m_smoothInput == SET_INDIVIDUALLY )
	{
		m_smoothcost  = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_labels];
		terminateOnError(!m_smoothcost,"Not enough memory");
		for ( i = 0; i < m_num_labels*m_num_labels; i++ ) 
			m_smoothcost[i] = (EnergyTermType) 0;

		m_smoothType      = ARRAY;
		m_varying_weights = 0;
	}
	else m_smoothType  = NONE;


	initialize_memory();
	//srand(time(NULL));
	srand(0);
    // <!-- bagon
    // sign class as valid
#ifdef MEX_COMPILE
    //class_sig = VALID_CLASS_SIGNITURE;
    class_sig = 0xabcd0123;
#else
    class_sig = 0xabcd0123;
#endif
    // bagon -->

	//<!-- Msun
	m_solver = GCsolver;
	m_nHigher = 0; // by default no higher-order-terms
	m_nHigherDet = 0; // by default no higher-order-terms
	m_nHigherOcc = 0; // by default no higher-order-terms
	m_nPairThing = 0; // by default no pair_thing terms
	m_useClassCoFlag = false; // by default ClassCo is not used
	m_num_comp_labels = 0;
    // 3/5 added
    m_seg_weight = 1;
	// Msun -->
}

/**************************************************************************************/
/* Use this constructor only for grid graphs                                          */
GCoptimization::GCoptimization(PixelType width,PixelType height,int nLabels,int dataSetup, int smoothSetup )
{
	commonGridInitialization(width,height,nLabels);
	
	m_labeling           = (LabelType *) new LabelType[m_num_pixels];
	terminateOnError(!m_labeling,"out of memory");
	for ( int i = 0; i < m_num_pixels; i++ ) m_labeling[i] = (LabelType) 0;

	// Msun: 2/5
	m_classExistFlag = new bool[m_num_labels];
	for ( int i = 0; i < m_num_labels; i++ ) m_classExistFlag[i] = false;
	m_classExistFlag[0] = true;
	// END	

	m_deleteLabeling = 1;
	commonInitialization(dataSetup,smoothSetup);	
}

/**************************************************************************************/
/* Use this constructor only for grid graphs                                          */
GCoptimization::GCoptimization(LabelType *m_answer,PixelType width,PixelType height,int nLabels,
							   int dataSetup, int smoothSetup) 
						   
{
	commonGridInitialization(width,height,nLabels);

	m_labeling = m_answer;
	

	for ( int i = 0; i < m_num_pixels; i++ )
		terminateOnError(m_labeling[i] < 0 || m_labeling[i] >= nLabels,"Initial labels are out of valid range");
	
	// Msun: 2/5
	m_classExistFlag = new bool[m_num_labels];
	for ( int i = 0; i < m_num_labels; i++ ) m_classExistFlag[i] = false;
	for ( int i = 0; i < m_num_pixels; i++ ){
		m_classExistFlag[ m_labeling[i]] = true;
	}
	// END	

	m_deleteLabeling = 0;
		
	commonInitialization(dataSetup,smoothSetup);	
}

/**************************************************************************************/
/* Use this constructor for general graphs                                            */
GCoptimization::GCoptimization(PixelType nupixels,int nLabels,int dataSetup, int smoothSetup )
{

	commonNonGridInitialization(nupixels, nLabels);
	
	m_labeling           = (LabelType *) new LabelType[m_num_pixels];
	terminateOnError(!m_labeling,"out of memory");
	for ( int i = 0; i < nupixels; i++ ) m_labeling[i] = (LabelType) 0;

	// Byung - debug
	if(nupixels != m_num_pixels) {
		printf("ERROR. Something wrong at nupixels != m_num_pixels\n");
		exit(1);
	}


	// Msun: 2/5
	m_classExistFlag = new bool[m_num_labels];
	for ( int i = 0; i < m_num_labels; i++ ) m_classExistFlag[i] = false;
	m_classExistFlag[0] = true;
	// END	

	m_deleteLabeling = 1;

	commonInitialization(dataSetup,smoothSetup);	

#ifdef LOG_ON
	std::cout<<"m_num_pixels="<<m_num_pixels<<std::endl;
#endif
}

/**************************************************************************************/
/* Use this constructor for general graphs                                            */
GCoptimization::GCoptimization(LabelType *m_answer, PixelType nupixels,int nLabels,int dataSetup, int smoothSetup)
{
	commonNonGridInitialization(nupixels, nLabels);
	

	m_labeling = m_answer;
	for ( int i = 0; i < m_num_pixels; i++ )
		terminateOnError(m_labeling[i] < 0 || m_labeling[i] >= nLabels,"Initial labels are out of valid range");

	// Msun: 2/5
	m_classExistFlag = new bool[m_num_labels];
	for ( int i = 0; i < m_num_labels; i++ ) m_classExistFlag[i] = false;
	for ( int i = 0; i < m_num_pixels; i++ ){
		m_classExistFlag[ m_labeling[i]] = true;
	}
	// END	

	m_deleteLabeling = 0;

	commonInitialization(dataSetup,smoothSetup);	
}

// <!--Msun initialize the CRF with higher order terms
/**************************************************************************************/
// general HOP
void GCoptimization::InitHOP(int nhigher)
{
	// set up higher-order terms
	m_nHigher = nhigher;
	//mexPrintf("m_nHigher=%d\n",m_nHigher);
	m_higherCost = new EnergyTermType[nhigher * (m_num_labels + 1)]; 
	m_higherElements = new int[nhigher];                                                         
	memset(m_higherElements, 0, nhigher * sizeof(int)); // make sure they are zero               
	m_higherTruncation = new EnergyTermType[nhigher];                                                  
	m_higherIndex = new int *[nhigher];                                                          
	memset(m_higherIndex, 0, nhigher * sizeof(int *));                                           
	m_higherWeights = new EnergyTermType *[nhigher]; // BAGON                                          
	memset(m_higherWeights, 0, nhigher * sizeof(EnergyTermType*)); // BAGON                            
	m_higherP = new EnergyTermType[nhigher];
}

// detection HOP
void GCoptimization::InitDetHOP(int nhigher)
{
	// set up higher-order terms
	m_nHigherDet = nhigher;
	
	// New concept -------------
	m_det_bool_label = new bool[nhigher];
	m_det_bool_negation = new bool[nhigher];
	memset( m_det_bool_label, true, nhigher * sizeof(bool)); // make sure they are true
	//if (!m_det_bool_label[0])
    //	mexPrintf("m_det_bool_label[0]=false\n");
	memset( m_det_bool_negation, false, nhigher * sizeof(bool)); // make sure they are zero
	//if (m_det_bool_negation[0])
    //	mexPrintf("m_det_bool_negation[0]=true\n");
	m_det_bool_lookupValid = (PixelType *) new PixelType[nhigher];
	m_det_fix_bool = new bool[nhigher];
	for (int i=0; i<nhigher; i++)
	{
		m_det_bool_lookupValid[i] = -1;
		m_det_fix_bool[i] = false;
	}
	
 	// -------------------------
	//mexPrintf("m_nHigher=%d\n",m_nHigher);
	m_higherDetCost = new EnergyTermType[nhigher * (m_num_labels + 1)];
	m_higherDetUW = new EnergyTermType[nhigher];
	m_det_labeling = new LabelType[nhigher];

	m_higherDetElements = new int[nhigher];
	memset(m_higherDetElements, 0, nhigher * sizeof(int)); // make sure they are zero               
	m_higherDetTruncation = new EnergyTermType[nhigher];                                                  
	m_higherDetIndex = new int *[nhigher]; 
	memset(m_higherDetIndex, 0, nhigher * sizeof(int *));                                           
	m_higherDetWeights = new EnergyTermType *[nhigher];
	memset(m_higherDetWeights, 0, nhigher * sizeof(EnergyTermType*));

	m_higherDetP = new EnergyTermType[nhigher];

    // 3/5 added class-specific weights
    m_det_class_weights = new EnergyTermType [m_num_labels];
	memset( m_det_class_weights, 0, m_num_labels * sizeof(EnergyTermType));
	for (int ii=0; ii<m_num_labels; ii++)
		m_det_class_weights[ii] = 1.0;
}

// occlusion HOP
void GCoptimization::InitOccHOP(int nhigher)
{
	// set up higher-order terms
	m_nHigherOcc = nhigher;
	//mexPrintf("m_nHigher=%d\n",m_nHigher);
	// New concept ----

	// <--Byung
	//m_det_bool_label = new bool[nhigher];
	//memset( m_det_bool_label, true, nhigher * sizeof(bool)); // make sure they are zero
	// Byung -->

	// <--Byung
	m_Occl2HigherDetIndex = new int[nhigher]; 
	memset(m_Occl2HigherDetIndex, 0, nhigher * sizeof(int)); // make sure they are zero
	// Byung -->


    // ----------------
	m_higherOccCost = new EnergyTermType[nhigher * (m_num_labels + 1)]; 
	m_higherOccElements = new int[nhigher];                                                         
	memset(m_higherOccElements, 0, nhigher * sizeof(int)); // make sure they are zero               
	m_higherOccTruncation = new EnergyTermType[nhigher];                                                  
	m_higherOccIndex = new int *[nhigher];                                                          
	memset(m_higherOccIndex, 0, nhigher * sizeof(int *));                                           
	m_higherOccWeights = new EnergyTermType *[nhigher]; // BAGON                                          
	memset(m_higherOccWeights, 0, nhigher * sizeof(EnergyTermType*)); // BAGON                            
	m_higherOccP = new EnergyTermType[nhigher];
}

// Pair_thing
void GCoptimization::InitPairThing(int nPairThing)
{
	m_nPairThing = nPairThing;

	m_PairThingIndice = new int *[nPairThing]; 
	memset( m_PairThingIndice, 0, nPairThing * sizeof(int *));

	m_PairThingWeights = new EnergyTermType *[nPairThing];
	memset( m_PairThingWeights, 0, nPairThing * sizeof(EnergyTermType*));

	m_random_pair_thing_order = true;
	m_PairThingTable = new int[nPairThing];
}

// class co-occurrence 
void GCoptimization::InitClassCo(bool useClassCoFlag)
{
	m_useClassCoFlag = useClassCoFlag;
//	if (m_useClassCoFlag){
//	}
}

// Msun -->

/**************************************************************************************/

void GCoptimization::setData(EnergyTermType* dataArray)
{
	terminateOnError(m_dataType != NONE,
		            "ERROR: you already set the data, or said you'll use member function setDataCost() to set data");	

    // <!-- bagon
    /* allocate memory dinamiocally */
    m_datacost = (EnergyTermType *) new EnergyTermType[m_num_pixels*m_num_labels];
    for ( int i(0); i < m_num_pixels*m_num_labels; i++ ){
        m_datacost[i] = dataArray[i];
	}
	// m_datacost  = dataArray;
    // bagon -->
	m_dataType  = ARRAY;
}

/**************************************************************************************/

void GCoptimization::setData(dataFnPix dataFn)
{
	terminateOnError(m_dataType != NONE,
		            "ERROR: you already set the data, or said you'll use member function setDataCost() to set data");	
	
	m_dataFnPix = dataFn;
	m_dataType  = FUNCTION_PIX;
}

/**************************************************************************************/

void GCoptimization::setData(dataFnCoord dataFn)
{
	terminateOnError(m_dataType != NONE,
		            "ERROR: you already set the data, or said you'll use member function setDataCost() to set data");	

	terminateOnError( !m_grid_graph,"Cannot use data function based on coordinates for non-grid graph");

	m_dataFnCoord = dataFn;
	m_dataType    = FUNCTION_COORD;
}

/**************************************************************************************/

void GCoptimization::setSmoothness(EnergyTermType* V)
{


#ifdef LOG_ON
	printf("We called correct setSmoothness\n"); fflush(stdout);
#endif
	terminateOnError(m_smoothType != NONE,
		            "ERROR: you already set smoothness, or said you'll use member function setSmoothCost() to set Smoothness Costs");	


	m_smoothType = ARRAY;
    // <!-- bagon
    /* allocate memory dinamiocally */
    m_smoothcost = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_labels];
    for ( int i(0); i < m_num_labels*m_num_labels; i++ )
        m_smoothcost[i] = V[i];
	// m_smoothcost = V;
    // bagon -->
    m_varying_weights = 0;

}

/**************************************************************************************/
void GCoptimization::setSmoothness(EnergyTermType* V,EnergyTermType* hCue, EnergyTermType* vCue)
{
	terminateOnError(m_smoothType != NONE,
		            "ERROR: you already set smoothness, or said you'll use member function setSmoothCost() to set Smoothness Costs");	

	terminateOnError(!m_grid_graph,
		            "ERROR: for a grid graph, you can't use vertical and horizontal cues.  Use setNeighbors() member function to encode spatially varying cues");	

	m_varying_weights = 1;

    m_smoothType      = ARRAY;
    // <!-- bagon
    /* allocate memory dinamiocally */
    int i = 0;
    m_smoothcost = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_labels];
    for ( i = 0; i < m_num_labels*m_num_labels; i++ )
        m_smoothcost[i] = V[i];
    
    m_vertWeights = (EnergyTermType *) new EnergyTermType[m_num_pixels];
    for ( i = 0; i < m_num_pixels; i++ )
        m_vertWeights[i] = vCue[i];
    
    m_horizWeights = (EnergyTermType *) new EnergyTermType[m_num_pixels];
    for ( i = 0; i < m_num_pixels; i++ )
        m_horizWeights[i] = hCue[i];
    
	//m_vertWeights     = vCue;
	//m_horizWeights    = hCue;
	//m_smoothcost      = V;
    // bagon -->
}

/**************************************************************************************/

void GCoptimization::setSmoothness(smoothFnCoord horz_cost, smoothFnCoord vert_cost)
{

	terminateOnError(m_smoothType != NONE,
		            "ERROR: you already set smoothness, or said you'll use member function setSmoothCost() to set Smoothness Costs");	

	terminateOnError( !m_grid_graph,"Cannot use smoothness function based on coordinates for non-grid graph");

	m_smoothType    = FUNCTION_COORD;
	m_horz_cost     = horz_cost;
	m_vert_cost     = vert_cost;
}

/**************************************************************************************/

void GCoptimization::setSmoothness(smoothFnPix cost) // Msun: this is more flexible... need to use this in the future
{
	terminateOnError(m_smoothType != NONE,
		            "ERROR: you already set smoothness, or said you'll use member function setSmoothCost() to set Smoothness Costs");	

	m_smoothType    = FUNCTION_PIX;
	m_smoothFnPix   = cost;
}

/**************************************************************************************/
// <!--Msun set the higher-order potential model terms
// set SegLoss
void GCoptimization::setSegLoss(EnergyTermType* lossArray)
{
    /* allocate memory dinamiocally */
    m_losscost = (EnergyTermType *) new EnergyTermType[m_num_pixels*m_num_labels];
    for ( int i(0); i < m_num_pixels*m_num_labels; i++ )
        m_losscost[i] = lossArray[i];
}

// printThing
void GCoptimization::printThing()
{
	for (int i=0; i<m_nHigherDet; i++)
		if (m_det_bool_label[i] == true){
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("1 ");
#endif
#endif
		}else{
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("0 ");
#endif
#endif
		}
}

// for general HOP terms
int GCoptimization::setOneHOP(int n, // number of nodes participating in this potential
            int* ind, // indices of participating nodes - needs to be allocated outside
            EnergyTermType* weights,  // w_i for each node in the potential - needs to be allocated outside
            EnergyTermType* gammas,  // gamma_l for all labels and gamma_max - do not allocate
            EnergyTermType Q) 
{
   	// find a vacant slot to insert this HOpotential
	int hoi(0), ii(0);
	for ( ; hoi < m_nHigher && m_higherElements[hoi] > 0 ; hoi++);
	if ( hoi == m_nHigher )
		return -1; // no vacant slot for this HOpotential

	m_higherElements[hoi] = n;
	m_higherTruncation[hoi] = Q;
	m_higherP[hoi] = 0;

	// do allocate

	m_higherWeights[hoi] = new EnergyTermType[n];
	m_higherIndex[hoi] = new int[n];

	for (ii = 0; ii<n ; ii++) {
		m_higherP[hoi]+=weights[ii];
		m_higherWeights[hoi][ii] = weights[ii];

		m_higherIndex[hoi][ii] = ind[ii];
	}
	if (2*Q >= m_higherP[hoi]) {
		// 2*Q must be smaller than P (see sec. 4.2 in tech-report)
		printf("2Q:%f,sum_w:%f\n", 2*Q, m_higherP[hoi]); fflush(stdout);
		return -1;
	}
	memcpy(m_higherCost + (m_num_labels+1)*hoi, gammas, (m_num_labels+1)*sizeof(EnergyTermType));

	return 1;
}

// for Detection HOP
int GCoptimization::setOneDetHOP(int n, // number of nodes participating in this potential
            int* ind, // indices of participating nodes - needs to be allocated outside
            EnergyTermType* weights,  // w_i for each node in the potential - needs to be allocated outside
            EnergyTermType* gammas,  // gamma_l for all labels and gamma_max - do not allocate
            EnergyTermType Q,
            EnergyTermType thing_uw)
{
   	// find a vacant slot to insert this HOpotential
	int hoi(0), ii(0);
	for ( ; hoi < m_nHigherDet && m_higherDetElements[hoi] > 0 ; hoi++);
	if ( hoi == m_nHigherDet )
		return -1; // no vacant slot for this HOpotential

	m_higherDetElements[hoi] = n;
	m_higherDetTruncation[hoi] = Q;
	m_higherDetP[hoi] = 0;
	m_higherDetUW[hoi] = thing_uw;
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("thing_uw[%d]=%f ",hoi,m_higherDetUW[hoi]);
#endif
#endif

	// do allocate

	m_higherDetWeights[hoi] = new EnergyTermType[n];
	m_higherDetIndex[hoi] = new int[n];

	for (ii = 0; ii<n ; ii++) {
		m_higherDetP[hoi]+=weights[ii];
		m_higherDetWeights[hoi][ii] = weights[ii];

		m_higherDetIndex[hoi][ii] = ind[ii];
	}
	//if (2*Q >= m_higherDetP[hoi]) {
		// 2*Q must be smaller than P (see sec. 4.2 in tech-report)
	//	printf("2Q:%f,sum_w:%f\n", 2*Q, m_higherP[hoi]); fflush(stdout);
	//	return -1;
	//}
	memcpy(m_higherDetCost + (m_num_labels+1)*hoi, gammas, (m_num_labels+1)*sizeof(EnergyTermType));

	// 2/2 find the first zeros gamma, which will be the det_label
	for (ii = 0; ii<m_num_labels ; ii++)
		if (gammas[ii] == 0){
			m_det_labeling[hoi] = ii;
			break;
		}
	

	return 1;
}

// for Occlusion HOP
int GCoptimization::setOneOccHOP(int n, // number of nodes participating in this potential
            int* ind, // indices of participating nodes - needs to be allocated outside
            EnergyTermType* weights,  // w_i for each node in the potential - needs to be allocated outside
            EnergyTermType* gammas,  // gamma_l for all labels and gamma_max - do not allocate
            EnergyTermType Q,
            int ind2Det) // indices of participating nodes - needs to be allocated outside
{
   	// find a vacant slot to insert this HOpotential
	int hoi(0), ii(0);
	for ( ; hoi < m_nHigherOcc && m_higherOccElements[hoi] > 0 ; hoi++);
	if ( hoi == m_nHigherOcc )
		return -1; // no vacant slot for this HOpotential

	// New concept ===
	m_Occl2HigherDetIndex[hoi] = ind2Det;
    // ===============

	m_higherOccElements[hoi] = n;
	m_higherOccTruncation[hoi] = Q;
	m_higherOccP[hoi] = 0;

	// do allocate

	m_higherOccWeights[hoi] = new EnergyTermType[n];
	m_higherOccIndex[hoi] = new int[n];

	for (ii = 0; ii<n ; ii++) {
		m_higherOccP[hoi]+=weights[ii];
		m_higherOccWeights[hoi][ii] = weights[ii];

		m_higherOccIndex[hoi][ii] = ind[ii];
	}
	if (2*Q >= m_higherOccP[hoi]) {
		// 2*Q must be smaller than P (see sec. 4.2 in tech-report)
		printf("2Q:%f,sum_w:%f\n", 2*Q, m_higherP[hoi]); fflush(stdout);
		return -1;
	}
	memcpy(m_higherOccCost + (m_num_labels+1)*hoi, gammas, (m_num_labels+1)*sizeof(EnergyTermType));

	return 1;
}

// pair-thing potential
void GCoptimization::setOnePairThing( int pti,
			int* ind,
			EnergyTermType* weights)
{
	m_PairThingIndice[pti] = new int[2];
	m_PairThingWeights[pti] = new EnergyTermType[4];

	m_PairThingIndice[pti][0] = ind[0];
	m_PairThingIndice[pti][1] = ind[1];

	m_PairThingWeights[pti][0] = weights[0];
	m_PairThingWeights[pti][1] = weights[1];
	m_PairThingWeights[pti][2] = weights[2];
	m_PairThingWeights[pti][3] = weights[3];

	for ( int ii=0; ii< m_nPairThing; ii++)
		m_PairThingTable[ii] = ii;
}

// class Co-occurrence
void GCoptimization::setClassCo(  EnergyTermType * uw, EnergyTermType * pw)
{
	if (!m_useClassCoFlag)
		return;

	int i;
	m_class_co_u_w = (EnergyTermType *) new EnergyTermType[m_num_labels];
    for ( i = 0; i < m_num_labels; i++ ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("uw[%d]=%f ",i,uw[i]);
#endif
#endif
        m_class_co_u_w[i] = uw[i];
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("\n");
#endif
#endif

	m_class_co_p_w = (EnergyTermType *) new EnergyTermType[m_num_labels*m_num_labels];
    for ( i = 0; i < m_num_labels*m_num_labels; i++ ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("pw[%d]=%f ",i,pw[i]);
#endif
#endif
        m_class_co_p_w[i] = pw[i];
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("\n");
#endif
#endif
}

void GCoptimization::setCompClassCo(  EnergyTermType * uw, EnergyTermType * pw)
{
	if (!m_useClassCoFlag)
		return;

	int i;
	m_class_co_u_w = (EnergyTermType *) new EnergyTermType[m_num_comp_labels];
    for ( i = 0; i < m_num_comp_labels; i++ ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("uw[%d]=%f ",i,uw[i]);
#endif
#endif
        m_class_co_u_w[i] = uw[i];
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("\n");
#endif
#endif

	m_class_co_p_w = (EnergyTermType *) new EnergyTermType[m_num_comp_labels*m_num_comp_labels];
    for ( i = 0; i < m_num_comp_labels*m_num_comp_labels; i++ ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("pw[%d]=%f ",i,pw[i]);
#endif
#endif
        m_class_co_p_w[i] = pw[i];
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("\n");
#endif
#endif
}


// for CompLabel
// Initization 
void GCoptimization::InitCompLabel(int num_comp_labels)
{
	// set up higher-order terms
	m_num_comp_labels = num_comp_labels;
	m_num_dup_labels = new int[num_comp_labels];                                                         
	memset(m_num_dup_labels, 0, num_comp_labels * sizeof(int)); // make sure they are zero               
	m_comp_label2label = new LabelType *[num_comp_labels];                                                          
	memset(m_comp_label2label, 0, num_comp_labels * sizeof(LabelType *));                                           
	m_label2comp_label = new LabelType [m_num_labels];                                                          
	memset(m_label2comp_label, 0, m_num_labels * sizeof(LabelType));                                           
}

// for comp_label2label matching
int GCoptimization::setOneCompLabel(int n, // number of nodes participating in this potential
            int* ind // indices of participating nodes - needs to be allocated outside
            ) 
{
   	// find a vacant slot to insert this HOpotential
	int compi(0), ii(0);
	for ( ; compi < m_num_comp_labels && m_num_dup_labels[compi] > 0 ; compi++);
	if ( compi == m_num_comp_labels )
		return -1; // no vacant slot for this HOpotential

	m_num_dup_labels[compi] = n;

	// do allocate
	m_comp_label2label[compi] = new int[n];

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("in setOneCompLabel");
#endif
#endif
	for (ii = 0; ii<n ; ii++) {
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf(" ind[%d]=%d ",ii,ind[ii]);
#endif
#endif
		m_comp_label2label[compi][ii] = ind[ii];
		m_label2comp_label[ind[ii]] = compi;
	}

	return 1;
}
/**************************************************************************************/
// for general HOP
GCoptimization::EnergyType GCoptimization::giveHOPEnergy()
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType *W = new GCoptimization::EnergyType[m_num_labels];
	GCoptimization::EnergyType he=0;

   	// collect HOpotenatials terms
	int i, j;
	for(i = 0; i < m_nHigher; i++)    // for each HOpotential
	{
		for(j = 0; j < m_num_labels; j++) W[j] = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherElements[i]; j++)
			W[m_labeling[ m_higherIndex[i][j]]]+=m_higherWeights[i][j];

		GCoptimization::EnergyType cost, minCost = m_higherCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max

		for(j = 0;j < m_num_labels; j++)
		{
			if( m_higherTruncation[i] ==0) exit(1);
			cost = m_higherCost[(m_num_labels + 1) * i + j] + // gamma_j 
				(m_higherP[i] - W[j])     //  P - sum w_i \delta_j(x_c)
					* (m_higherCost[(m_num_labels + 1) * i + m_num_labels]-m_higherCost[(m_num_labels + 1) * i + j]) // gamma_max - gamma_j
					* (1 / m_higherTruncation[i]);    // 1 / Q
			if (minCost >= cost) minCost = cost;
		}
		// add HOpotential's energy to the total term
		
		he += minCost;
	}

    // add seg_weight
    he *= m_seg_weight;
	delete [] W;
	return he;
}

// for class specific detecion HOP2 // for byung
GCoptimization::EnergyType GCoptimization::giveClassDetHOP2Energy( LabelType class_label )
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType W, he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		if ( m_det_labeling[i] != class_label)
			continue;

		W = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherDetElements[i]; j++){
			if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i]){
				W+=m_higherDetWeights[i][j];
			}
		}

		GCoptimization::EnergyType minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max
		GCoptimization::EnergyType thing_uw = m_higherDetUW[i];

		if (m_det_bool_label[i]){
			minCost = (m_higherDetP[i] - W)     //  P - sum w_i \delta_j(x_c)
					* gamma_max // gamma_max
					* (1 / m_higherDetTruncation[i]);    // 1 / Q
		}else{
			// if det indicator is flase, force to use gamma_max
#ifdef NO_UNARY_DET_COMPILE
			minCost = gamma_max*(1+(W/(m_higherDetP[i]-m_higherDetTruncation[i])));
#else
			minCost = gamma_max*((W/(m_higherDetP[i]-m_higherDetTruncation[i])))+thing_uw;
			//minCost = gamma_max*((W/(m_higherDetP[i]-m_higherDetTruncation[i])));
			//minCost = thing_uw;
#endif
		}
		// add HOpotential's energy to the total term
		
		he_det += minCost;
	}

	// add weights
    he_det *= m_det_class_weights[class_label];

	return he_det;
}

// only get potential related to the u(y,X)
GCoptimization::EnergyType GCoptimization::giveClassDetHOP2XYEnergy( LabelType class_label )
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType W, he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		if ( m_det_labeling[i] != class_label)
			continue;

		W = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherDetElements[i]; j++){
			if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i]){
				W+=m_higherDetWeights[i][j];
			}
		}

		GCoptimization::EnergyType minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max
		GCoptimization::EnergyType thing_uw = m_higherDetUW[i];

		if (m_det_bool_label[i]){
			minCost = (m_higherDetP[i] - W)     //  P - sum w_i \delta_j(x_c)
					* gamma_max // gamma_max
					* (1 / m_higherDetTruncation[i]);    // 1 / Q
		}else{
			// if det indicator is flase, force to use gamma_max
			minCost = gamma_max*(W/(m_higherDetP[i]-m_higherDetTruncation[i]));
		}
		// add HOpotential's energy to the total term
		
		he_det += minCost;
	}

	// add weights
    he_det *= m_det_class_weights[class_label];

	return he_det;
}

// for class specific detecion HOP2 // for byung
GCoptimization::EnergyType GCoptimization::giveClassDetHOP2PriorEnergy( LabelType class_label )
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType W, he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		if ( m_det_labeling[i] != class_label)
			continue;

		GCoptimization::EnergyType minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max
		GCoptimization::EnergyType thing_uw = m_higherDetUW[i];

		if (!m_det_bool_label[i]){
#ifdef NO_UNARY_DET_COMPILE
			minCost = gamma_max;
#else
			minCost = thing_uw;
#endif
		}else{
			minCost = 0.0;
		}
		// add HOpotential's energy to the total term
		he_det += minCost;
	}

	// add weights
    he_det *= m_det_class_weights[class_label];

	return he_det;
}

// for detecion HOP2 // for byung
GCoptimization::EnergyType GCoptimization::giveDetHOP2Energy( )
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType W, he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		W = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherDetElements[i]; j++)
			if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i])
				W+=m_higherDetWeights[i][j];

		GCoptimization::EnergyType minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max
		GCoptimization::EnergyType thing_uw = m_higherDetUW[i];

		if (m_det_bool_label[i]){
			minCost = (m_higherDetP[i] - W)     //  P - sum w_i \delta_j(x_c)
					* gamma_max // gamma_max
					* (1 / m_higherDetTruncation[i]);    // 1 / Q
			if (m_higherDetTruncation[i]==0) {
#ifdef LOG_ON
				printf("minCost:%f\n",minCost); fflush(stdout);
#endif
			}
		}else{
			// if det indicator is flase, force to use gamma_max
#ifdef NO_UNARY_DET_COMPILE
			minCost = gamma_max*(1+(W/(m_higherDetP[i]-m_higherDetTruncation[i])));
#else
			minCost = gamma_max*((W/(m_higherDetP[i]-m_higherDetTruncation[i])))+thing_uw;
#endif
			if(m_higherDetP[i]-m_higherDetTruncation[i] == 0) {
#ifdef LOG_ON
				printf("minCost:%f\n",minCost); fflush(stdout);
#endif
			}
		}
	
		// add weights
		minCost *= m_det_class_weights[m_det_labeling[i]];
#ifdef MEX_COMPILE
		//mexPrintf("DetClassWei[%d]=%f\n",m_det_labeling[i],m_det_class_weights[m_det_labeling[i]]);
#else
		//Printf("DetClassWei[%d]=%f\n",m_det_labeling[i],m_det_class_weights[m_det_labeling[i]]);
#endif
	
		// add HOpotential's energy to the total term
		he_det += minCost;
	}

	if (he_det == !he_det){
#ifdef LOG_ON
		printf("he_det:%f\n",he_det); fflush(stdout);
#endif
	}
	return he_det;
}

// for class specific detecion HOP // for byung
GCoptimization::EnergyType GCoptimization::giveClassDetHOPEnergy( LabelType class_label )
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType *W = new GCoptimization::EnergyType[m_num_labels];
	GCoptimization::EnergyType he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		if ( m_det_labeling[i] != class_label)
			continue;

		for(j = 0; j < m_num_labels; j++) W[j] = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherDetElements[i]; j++)
			W[m_labeling[ m_higherDetIndex[i][j]]]+=m_higherDetWeights[i][j];

		GCoptimization::EnergyType cost, minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max

		if (m_det_bool_label[i]){
			minCost = std::numeric_limits<GCoptimization::EnergyType>::max( );
			for(j = 0;j < m_num_labels; j++)
			{
				cost = m_higherDetCost[(m_num_labels + 1) * i + j] + // gamma_j 
					(m_higherDetP[i] - W[j])     //  P - sum w_i \delta_j(x_c)
						* (m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]-m_higherDetCost[(m_num_labels + 1) * i + j]) // gamma_max - gamma_j
						* (1 / m_higherDetTruncation[i]);    // 1 / Q
				if (minCost >= cost) minCost = cost;
			}
		}else{
			// if det indicator is flase, force to use gamma_max
			minCost = gamma_max;
		}
		
		// add HOpotential's energy to the total term
		he_det += minCost;
	}

	delete [] W;

	// add weights
    he_det *= m_det_class_weights[class_label];

	return he_det;
}

// for detecion HOP
GCoptimization::EnergyType GCoptimization::giveDetHOPEnergy()
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType *W = new GCoptimization::EnergyType[m_num_labels];
	GCoptimization::EnergyType he_det=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherDet=%d", m_nHigherDet); fflush(stdout);
	for(i = 0; i < m_nHigherDet; i++)    // for each HOpotential
	{
		for(j = 0; j < m_num_labels; j++) W[j] = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherDetElements[i]; j++)
			W[m_labeling[ m_higherDetIndex[i][j]]]+=m_higherDetWeights[i][j];

		GCoptimization::EnergyType cost, minCost, gamma_max = m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max

		if (m_det_bool_label[i]){
			minCost = std::numeric_limits<GCoptimization::EnergyType>::max( );
			for(j = 0;j < m_num_labels; j++)
			{
				cost = m_higherDetCost[(m_num_labels + 1) * i + j] + // gamma_j 
					(m_higherDetP[i] - W[j])     //  P - sum w_i \delta_j(x_c)
						* (m_higherDetCost[(m_num_labels + 1) * i + m_num_labels]-m_higherDetCost[(m_num_labels + 1) * i + j]) // gamma_max - gamma_j
						* (1 / m_higherDetTruncation[i]);    // 1 / Q
				if (minCost >= cost) minCost = cost;
			}
		}else{
			// if det indicator is flase, force to use gamma_max
			minCost = gamma_max;
		}

		// add weights
		minCost *= m_det_class_weights[m_det_labeling[i]];
		
		// add HOpotential's energy to the total term
		he_det += minCost;
	}
	delete [] W;

	return he_det;
}

// for class specific occlusion HOP // for byung
GCoptimization::EnergyType GCoptimization::giveClassOccHOPEnergy( LabelType class_label)
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType *W = new GCoptimization::EnergyType[m_num_labels];
	GCoptimization::EnergyType he_occ=0;

   	// collect HOpotenatials terms
	int i, j;
	//printf("debug:m_nHigherOcc=%d", m_nHigherOcc); fflush(stdout);
	for(i = 0; i < m_nHigherOcc; i++)    // for each HOpotential
	{
		if ( m_det_labeling[ m_Occl2HigherDetIndex[i] ] != class_label)
			continue;

		for(j = 0; j < m_num_labels; j++) W[j] = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherOccElements[i]; j++)
			W[m_labeling[ m_higherOccIndex[i][j]]]+=m_higherOccWeights[i][j];

		GCoptimization::EnergyType cost, minCost, gamma_max = m_higherOccCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max

		if (m_det_bool_label[ m_Occl2HigherDetIndex[i] ]){
			minCost = std::numeric_limits<GCoptimization::EnergyType>::max( );
			for(j = 0;j < m_num_labels; j++)
			{
				cost = m_higherOccCost[(m_num_labels + 1) * i + j] + // gamma_j 
					(m_higherOccP[i] - W[j])     //  P - sum w_i \delta_j(x_c)
						* (m_higherOccCost[(m_num_labels + 1) * i + m_num_labels]-m_higherOccCost[(m_num_labels + 1) * i + j]) // gamma_max - gamma_j
						* (1 / m_higherOccTruncation[i]);    // 1 / Q
				if (minCost >= cost) minCost = cost;
			}
		}else{
			// if det indicator is flase, force to use gamma_max
			minCost = gamma_max;
		}
		// add HOpotential's energy to the total term
		
		he_occ += minCost;
	}

	delete [] W;

	return he_occ;
}
// for occlusion HOP
GCoptimization::EnergyType GCoptimization::giveOccHOPEnergy()
{
	// sum w_i delta_j(x_i)
	GCoptimization::EnergyType *W = new GCoptimization::EnergyType[m_num_labels];
	GCoptimization::EnergyType he_occ=0;

   	// collect HOpotenatials terms
	int i, j;
	for(i = 0; i < m_nHigherOcc; i++)    // for each HOpotential
	{
		for(j = 0; j < m_num_labels; j++) W[j] = 0;

		// count how many nodes are labeled L in the potential i
		for(j = 0; j < m_higherOccElements[i]; j++)
			W[m_labeling[ m_higherOccIndex[i][j]]]+=m_higherOccWeights[i][j];

		GCoptimization::EnergyType cost, minCost, gamma_max = m_higherOccCost[(m_num_labels + 1) * i + m_num_labels]; // gamma_max

		if (m_det_bool_label[ m_Occl2HigherDetIndex[i] ]){
			minCost = std::numeric_limits<GCoptimization::EnergyType>::max( );
			for(j = 0;j < m_num_labels; j++)
			{
				cost = m_higherOccCost[(m_num_labels + 1) * i + j] + // gamma_j 
					(m_higherOccP[i] - W[j])     //  P - sum w_i \delta_j(x_c)
						* (m_higherOccCost[(m_num_labels + 1) * i + m_num_labels]-m_higherOccCost[(m_num_labels + 1) * i + j]) // gamma_max - gamma_j
						* (1 / m_higherOccTruncation[i]);    // 1 / Q
				if (minCost >= cost) minCost = cost;
			}
		}else{
			// if det indicator is flase, force to use gamma_max
			minCost = gamma_max;
		}
		// add HOpotential's energy to the total term
		
		he_occ += minCost;
	}
	delete [] W;
	return he_occ;
}

// for each specific pair-thing potential // for byung. It is a bit in-efficient
GCoptimization::EnergyType GCoptimization::giveSpecPairThingEnergy(int thing1, int thing2)
{
	// search for which m_nPairThing
	GCoptimization::EnergyType pe = 0;
	for (int ii=0; ii< m_nPairThing; ii++){
		if ((m_PairThingIndice[ii][0]==thing1 && m_PairThingIndice[ii][1] == thing2) ||
			(m_PairThingIndice[ii][0]==thing2 && m_PairThingIndice[ii][1] == thing1)){
			if ( m_det_bool_label[m_PairThingIndice[ii][0]] == false &&
				 m_det_bool_label[m_PairThingIndice[ii][1]] == false){
				pe += m_PairThingWeights[ii][0];
			}else if ( m_det_bool_label[m_PairThingIndice[ii][0]] == false &&
					   m_det_bool_label[m_PairThingIndice[ii][1]] == true){
				pe += m_PairThingWeights[ii][1];
			}else if (m_det_bool_label[m_PairThingIndice[ii][0]] == true &&
						m_det_bool_label[m_PairThingIndice[ii][1]] == false){
				pe += m_PairThingWeights[ii][2];
			}else{
				pe += m_PairThingWeights[ii][3];
			}
			break;
		}
	}
	return pe;
}

// for pair-thing potential
GCoptimization::EnergyType GCoptimization::givePairThingEnergy()
{
	GCoptimization::EnergyType pe = 0;
//	mexPrintf("\nm_nPairThing=%d\n",m_nPairThing);
	for (int ii=0; ii< m_nPairThing; ii++){
		if ( m_det_bool_label[m_PairThingIndice[ii][0]] == false &&
			 m_det_bool_label[m_PairThingIndice[ii][1]] == false){
//			mexPrintf("ii=%d y[%d]=0,y[%d]=0,pe=%f\n",ii,m_PairThingIndice[ii][0],m_PairThingIndice[ii][1],m_PairThingWeights[ii][0]);
			pe += m_PairThingWeights[ii][0];
		}else if ( m_det_bool_label[m_PairThingIndice[ii][0]] == false &&
				   m_det_bool_label[m_PairThingIndice[ii][1]] == true){
//			mexPrintf("ii=%d y[%d]=0,y[%d]=1,pe=%f\n",ii,m_PairThingIndice[ii][0],m_PairThingIndice[ii][1],m_PairThingWeights[ii][1]);
			pe += m_PairThingWeights[ii][1];
		}else if (m_det_bool_label[m_PairThingIndice[ii][0]] == true &&
					m_det_bool_label[m_PairThingIndice[ii][1]] == false){
//			mexPrintf("ii=%d y[%d]=1,y[%d]=0,pe=%f\n",ii,m_PairThingIndice[ii][0],m_PairThingIndice[ii][1],m_PairThingWeights[ii][2]);
			pe += m_PairThingWeights[ii][2];
		}else{
//			mexPrintf("ii=%d y[%d]=1,y[%d]=1,pe=%f\n",ii,m_PairThingIndice[ii][0],m_PairThingIndice[ii][1],m_PairThingWeights[ii][3]);
			pe += m_PairThingWeights[ii][3];
		}
	}
	return pe;
}

// for Class Co-occurrence
GCoptimization::EnergyType GCoptimization::giveClassCoEnergy()
{
	if (!m_useClassCoFlag){
		return 0;
	}

	int ii,jj;
	GCoptimization::EnergyType coe = 0;

//<--Byung : debug co-occur
	for (ii = 0; ii < m_num_labels; ii++){
		if ( m_classExistFlag[ii] ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("l=%d ",ii);
#endif
#endif
		}
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("\n\n");
#endif
#endif
//-->Byung

	for (ii = 0; ii < m_num_labels; ii++){
		if ( m_classExistFlag[ii] ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("l=%d,coe=%f ",ii,m_class_co_u_w[ii]);
#endif
#endif
			coe += m_class_co_u_w[ii];
			for (jj = ii+1; jj < m_num_labels; jj++)
				if ( m_classExistFlag[jj] ){
#ifdef MEX_COMPILE
#ifdef LOG_ON
					mexPrintf("l1=%d,l2=%d,coe=%f ",ii,jj,m_class_co_p_w[ii+jj*m_num_labels]);
#endif
#endif
					coe += m_class_co_p_w[ii+jj*m_num_labels];
				}
		}
	}
	return coe;
}

// for Class Co-occurrence
GCoptimization::EnergyType GCoptimization::giveCompClassCoEnergy()
{
	if (!m_useClassCoFlag){
	#ifdef MEX_COMPILE
	#ifdef LOG_ON
		mexPrintf("m_useClassCoFlag=false\n");
#endif
	#endif
		return 0;
	}
	#ifdef MEX_COMPILE
	#ifdef LOG_ON
	mexPrintf("m_num_comp_labels=%d\n",m_num_comp_labels);
#endif
	#endif

	int ii,jj,sub_ii,sub_jj;
	GCoptimization::EnergyType coe = 0;
	
	for (ii = 0; ii < m_num_comp_labels; ii++){
		for (sub_ii=0; sub_ii< m_num_dup_labels[ii]; sub_ii++){
			LabelType target_label_ii = m_comp_label2label[ii][sub_ii];
			if ( m_classExistFlag[target_label_ii]){
	#ifdef MEX_COMPILE
	#ifdef LOG_ON
				mexPrintf("l=%d,coe=%f ",ii,m_class_co_u_w[ii]);
#endif
	#endif
				coe += m_class_co_u_w[ii];
				for (jj = ii+1; jj < m_num_comp_labels; jj++){
					for (sub_jj=0; sub_jj< m_num_dup_labels[jj]; sub_jj++){
						LabelType target_label_jj = m_comp_label2label[jj][sub_jj];
						if ( m_classExistFlag[ target_label_jj]){
		#ifdef MEX_COMPILE
		#ifdef LOG_ON
							//mexPrintf("l1=%d,l2=%d,coe=%f ",ii,jj,m_class_co_p_w[ii+jj*m_num_labels]);
							mexPrintf("l1=%d,l2=%d,coe=%f ",ii,jj,m_class_co_p_w[ii+jj*m_num_comp_labels]);
#endif
		#endif
							//coe += m_class_co_p_w[ii+jj*m_num_labels];
							coe += m_class_co_p_w[ii+jj*m_num_comp_labels];
							break;
						}
					}
				}
				break;
			}
		}
	}
	return coe;
}

// get indicator labels
void 
GCoptimization::ExportDetLabels(LabelType* labels)
{
    for( int i(0); i < m_nHigherDet; i++){
		if (m_det_bool_label[i]){
	        labels[i] = 1;
		}else{
	        labels[i] = 0;
		}
	}
}

// set detection indicator labels // for byung
void 
GCoptimization::SetDetBoolLabels(const bool* det_bool_labels)
{
    for( int i(0); i < m_nHigherDet; i++ ) {
		m_det_bool_label[i] = det_bool_labels[i];
		//if (m_det_bool_label[i])
		//	mexPrintf("m_det_bool_label[%d]=true ",i);
    }
}

// Msun: 9/13 added fix detection indicators
void GCoptimization::SetDetFixBool(const bool* det_fix_bool)
{
    for( int i(0); i < m_nHigherDet; i++ ) {
		m_det_fix_bool[i] = det_fix_bool[i];
	}
}

// Msun: 3/5 added loss energy
GCoptimization::EnergyType GCoptimization::giveSegLossEnergyArray()
{
	EnergyType eng = (EnergyType) 0;


	for ( int i = 0; i < m_num_pixels; i++ ) if (m_labeling[i] != m_num_labels+1)
		eng = eng + m_losscost(i,m_labeling[i]);

	return(eng);
}
// Msun -->

/**************************************************************************************/
// <!-- bagon
void 
GCoptimization::SetAllLabels(const LabelType* labels)
{
	//<-- Byung 5/6 Debugging
	for (int ii(0); ii < m_num_labels; ii++){
		m_classExistFlag[ii] = false;
	}
	m_classExistFlag[0] = true; // void always exists
	//--> Byung

    for( int i(0); i < m_num_pixels; i++ ) {
//		mexPrintf("label=%d\n",labels[i]);
        terminateOnError( (labels[i]<0) || (labels[i]>=m_num_labels), "Wrong label value");
        m_labeling[i] = labels[i];
		//Msun 2/5
		m_classExistFlag[ m_labeling[i]] = true;
		//END
    }
}
/**************************************************************************************/    
void 
GCoptimization::ExportLabels(LabelType* labels)
{
    for( int i(0); i < m_num_pixels; i++)
        labels[i] = m_labeling[i];
}
// bagon -->
/**************************************************************************************/    
GCoptimization::EnergyType GCoptimization::giveDataEnergy()
{
	if ( m_dataType == ARRAY) 
		return(giveDataEnergyArray()*m_seg_weight);
	else if ( m_dataType == FUNCTION_PIX )
		return(giveDataEnergyFnPix()*m_seg_weight);
	else if (m_dataType == FUNCTION_COORD ) return(giveDataEnergyFnCoord()*m_seg_weight);
	else terminateOnError(1,"Did not initialize the data costs yet");

	return(0);
}

/**************************************************************************************/
	
GCoptimization::EnergyType GCoptimization::giveDataEnergyArray()
{
	EnergyType eng = (EnergyType) 0;


	for ( int i = 0; i < m_num_pixels; i++ ) if (m_labeling[i] != m_num_labels+1)
		eng = eng + m_datacost(i,m_labeling[i]);

	return(eng);
}

/**************************************************************************************/
	
GCoptimization::EnergyType GCoptimization::giveDataEnergyFnPix()
{

	EnergyType eng = (EnergyType) 0;

	for ( int i = 0; i < m_num_pixels; i++ )
		eng = eng + m_dataFnPix(i,m_labeling[i]);

	return(eng);
}
/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveDataEnergyFnCoord()
{
	EnergyType eng = (EnergyType) 0;

	for ( int y = 0; y < m_height; y++ )
		for ( int x = 0; x < m_width; x++ )
			eng = eng + m_dataFnCoord(x,y,m_labeling[x+y*m_width]);

	return(eng);

}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy()
{

	if ( m_grid_graph )
	{
		if ( m_smoothType == ARRAY )
		{
			if (m_varying_weights) return(giveSmoothEnergy_G_ARRAY_VW());
			else return(giveSmoothEnergy_G_ARRAY()*m_seg_weight);
		}
		else if ( m_smoothType == FUNCTION_PIX ) return(giveSmoothEnergy_G_FnPix()*m_seg_weight);
		else if ( m_smoothType == FUNCTION_COORD ) return(giveSmoothEnergy_G_FnCoord()*m_seg_weight);
		else terminateOnError(1,"Did not initialize smoothness costs yet, can't compute smooth energy");
	}
	else
	{
		if ( m_smoothType == ARRAY ) return(giveSmoothEnergy_NG_ARRAY()*m_seg_weight);
		else if ( m_smoothType == FUNCTION_PIX ) return(giveSmoothEnergy_NG_FnPix()*m_seg_weight);
		else terminateOnError(1,"Did not initialize smoothness costs yet, can't compute smooth energy");
	}
	return(0);

}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_NG_FnPix()
{

	EnergyType eng = (EnergyType) 0;
	int i;
	Neighbor *temp; 

	for ( i = 0; i < m_num_pixels; i++ )
		if ( !m_neighbors[i].isEmpty() )
		{
			m_neighbors[i].setCursorFront();
			while ( m_neighbors[i].hasNext() )
			{
				temp = (Neighbor *) m_neighbors[i].next();
				if ( i < temp->to_node )
					eng = eng + m_smoothFnPix(i,temp->to_node, m_labeling[i],m_labeling[temp->to_node]);
			}
		}
		
	return(eng);
	
}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_NG_ARRAY()
{
	EnergyType eng = (EnergyType) 0;
	int i;
	Neighbor *temp; 

	for ( i = 0; i < m_num_pixels; i++ )
		if ( !m_neighbors[i].isEmpty() )
		{
			m_neighbors[i].setCursorFront();
			while ( m_neighbors[i].hasNext() )
			{
				temp = (Neighbor *) m_neighbors[i].next();

				if ( i < temp->to_node )
					eng = eng + m_smoothcost(m_labeling[i],m_labeling[temp->to_node])*(temp->weight);

			}
		}

	return(eng);
}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_ARRAY_VW()
{

	EnergyType eng = (EnergyType) 0;
	int x,y,pix;

	// printf("\nIn right place");
	
	for ( y = 0; y < m_height; y++ )
		for ( x = 1; x < m_width; x++ ) if ((m_labeling[x+y*m_width-1] != m_num_labels +1) &&(m_labeling[x+y*m_width] != m_num_labels +1) ) 
		{
			pix = x+y*m_width;
			eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1])*m_horizWeights[pix-1];
		}

	for ( y = 1; y < m_height; y++ )
		for ( x = 0; x < m_width; x++ )if ((m_labeling[x+y*m_width-m_width] != m_num_labels +1) &&(m_labeling[x+y*m_width] != m_num_labels +1) ) 
		{
			pix = x+y*m_width;
			eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width])*m_vertWeights[pix-m_width];
		}

	
	return(eng);
	
}
/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_ARRAY()
{

	EnergyType eng = (EnergyType) 0;
	int x,y,pix;


	for ( y = 0; y < m_height; y++ )
		for ( x = 1; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1]);
		}

	for ( y = 1; y < m_height; y++ )
		for ( x = 0; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width]);
		}

	return(eng);
	
}
/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_FnPix()
{

	EnergyType eng = (EnergyType) 0;
	int x,y,pix;


	for ( y = 0; y < m_height; y++ )
		for ( x = 1; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_smoothFnPix(pix,pix-1,m_labeling[pix],m_labeling[pix-1]);
		}

	for ( y = 1; y < m_height; y++ )
		for ( x = 0; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_smoothFnPix(pix,pix-m_width,m_labeling[pix],m_labeling[pix-m_width]);
		}

	return(eng);
	
}
/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_FnCoord()
{

	EnergyType eng = (EnergyType) 0;
	int x,y,pix;


	for ( y = 0; y < m_height; y++ )
		for ( x = 1; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_horz_cost(x-1,y,m_labeling[pix],m_labeling[pix-1]);
		}

	for ( y = 1; y < m_height; y++ )
		for ( x = 0; x < m_width; x++ )
		{
			pix = x+y*m_width;
			eng = eng + m_vert_cost(x,y-1,m_labeling[pix],m_labeling[pix-m_width]);
		}

	return(eng);
	
}

/**************************************************************************************/
 
GCoptimization::EnergyType GCoptimization::compute_energy()
{
	//printf("Data:%f,Smooth:%f,HOP:%f,DetHOP:%f,OccHOP:%f,PairTh:%f,ClassCo:%f\n", giveDataEnergy(),giveSmoothEnergy(),giveHOPEnergy(),giveDetHOPEnergy(),giveOccHOPEnergy(),givePairThingEnergy(),giveClassCoEnergy()); fflush(stdout);
	/*printf("in compute_energy()\n"); fflush(stdout);
	printf("before compute_data_energy()\n"); fflush(stdout);
    giveDataEnergy();
	printf("before compute_sn_energy()\n"); fflush(stdout);
	giveSmoothEnergy();
	printf("before compute_HOP_energy()\n"); fflush(stdout);
	giveHOPEnergy();
	printf("before compute_DetHOP_energy()\n"); fflush(stdout);
	giveDetHOPEnergy();*/
  	EnergyType class_co_en = 0;
	if (m_num_comp_labels == 0){
		#ifdef MEX_COMPILE
		#ifdef LOG_ON
		mexPrintf("m_num_comp_labels == 0");
#endif
		#endif
		class_co_en = giveClassCoEnergy();
	}else{
		#ifdef MEX_COMPILE
		#ifdef LOG_ON
		mexPrintf("m_num_comp_labels = %d\n",m_num_comp_labels);
#endif
		#endif
		class_co_en = giveCompClassCoEnergy();
	}

#ifdef OLD_DET_COMPILE
	return(giveDataEnergy()+giveSmoothEnergy()+
			giveHOPEnergy()+giveDetHOPEnergy()+giveOccHOPEnergy()+
			givePairThingEnergy()
			+class_co_en+
            giveSegLossEnergyArray());
#else
	return(giveDataEnergy()+giveSmoothEnergy()+
			giveHOPEnergy()+giveDetHOP2Energy()+giveOccHOPEnergy()+
			givePairThingEnergy()
			+class_co_en+
            giveSegLossEnergyArray());
#endif

}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::expansion(int max_num_iterations)
{
	return(start_expansion(max_num_iterations)); 
}

/**************************************************************************************/


GCoptimization::EnergyType GCoptimization::expansion()
{
	return(start_expansion(MAX_INTT));
}

/**************************************************************************************/


GCoptimization::EnergyType GCoptimization::start_expansion(int max_num_iterations )
{
	
	int curr_cycle = 1;
	EnergyType new_energy,old_energy;
	

	new_energy = compute_energy();

	//old_energy = (new_energy+1)*10; // BAGON changed init value to exceed current energy by factor of 10 (thanks to A. Khan)
	old_energy = (new_energy+1000); // Msun: now energy might be negtive due to the loss function

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
#ifdef LOG_ON
		printf("debug:before one expansion\n"); fflush(stdout);
#endif
		new_energy = oneExpansionIteration();
		
		curr_cycle++;	
	}

	return(new_energy);
}

/****************************************************************************/

void GCoptimization::scramble_label_table()
{

#ifdef LOG_ON
	printf("random alpha expansion!\n"); fflush(stdout);
#endif
   LabelType r1,r2,temp;
   int num_times,cnt;


   num_times = m_num_labels*2;

   for ( cnt = 0; cnt < num_times; cnt++ )
   {
      r1 = rand()%m_num_labels;  
      r2 = rand()%m_num_labels;  

      temp             = m_labelTable[r1];
      m_labelTable[r1] = m_labelTable[r2];
      m_labelTable[r2] = temp;
   }
}

//<--Msun

void GCoptimization::scramble_pair_thing_table()
{
   int r1,r2,temp;
   int num_times,cnt;
#ifdef MEX_COMPILE
#ifdef LOG_ON
   mexPrintf("in scramble_pair_thing_table");
#endif
#endif

#ifdef LOG_ON
	printf("random alpha expansion!\n"); fflush(stdout);
#endif

   num_times = m_nPairThing*2;

   for ( cnt = 0; cnt < num_times; cnt++ )
   {
      r1 = rand()%m_nPairThing;  
      r2 = rand()%m_nPairThing;  

      temp             = m_PairThingTable[r1];
      m_PairThingTable[r1] = m_PairThingTable[r2];
      m_PairThingTable[r2] = temp;
   }
}

//-->

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::alpha_expansion(LabelType label)
{
	terminateOnError( label < 0 || label >= m_num_labels,"Illegal Label to Expand On");
	
	if (m_solver == GCsolver) perform_alpha_expansion(label);
	else if (m_solver == QPBOsolver) perform_alpha_expansion_QPBO(label);
	else terminateOnError(1,"Did not specify the solver!");
    
	return(compute_energy());
}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::alpha_expansion(LabelType alpha_label, PixelType *pixels, int num )
{
	PixelType i,size = 0; 
#ifdef MEX_COMPILE
	Energy *e = new Energy(mexErrMsgTxt);
#else
	Energy *e = new Energy();
#endif


	for ( i = 0; i<num ; i++ )
	{
		if ( m_labeling[pixels[i]] != alpha_label )
		{
			m_lookupPixVar[pixels[i]] = i;
		}
	}

	
	if ( size > 0 ) 
	{
		Energy::Var *variables = (Energy::Var *) new Energy::Var[size];

		for ( i = 0; i < size; i++ )
			variables[i] = e ->add_variable();

		if ( m_dataType == ARRAY ) add_t_links_ARRAY(e,variables,size,alpha_label);
		else  if  ( m_dataType == FUNCTION_PIX ) add_t_links_FnPix(e,variables,size,alpha_label);
		else  add_t_links_FnCoord(e,variables,size,alpha_label);


		if ( m_grid_graph )
		{
			if ( m_smoothType == ARRAY )
			{
				if (m_varying_weights) set_up_expansion_energy_G_ARRAY_VW_pix(size,alpha_label,e,variables,pixels,num);
				else set_up_expansion_energy_G_ARRAY_pix(size,alpha_label,e,variables,pixels,num);
			}
			else if ( m_smoothType == FUNCTION_PIX ) {
#ifdef MEX_COMPILE
                mexErrMsgIdAndTxt("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
#else
                printf("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
				exit(1);
#endif
//				printf("NOT SUPPORTED YET,exiting!");
//				exit(0);
				//if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_G_FnPix_pix(size,alpha_label,e,variables);
			}
			else
			{
#ifdef MEX_COMPILE
                mexErrMsgIdAndTxt("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
#else
                printf("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
				exit(1);
#endif
//				printf("NOT SUPPORTED YET,exiting!");
//				exit(0);
				//set_up_expansion_energy_G_FnCoord_pix(size,alpha_label,e,variables);
			}
			
		}
		else
		{
#ifdef MEX_COMPILE
			mexErrMsgIdAndTxt("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
#else
			printf("GraphCut:internal", "NOT SUPPORTED YET,exiting!");
			exit(1);
#endif

//			printf("NOT SUPPORTED YET,exiting!");
//			exit(0);
			/*if ( m_smoothType == ARRAY ) set_up_expansion_energy_NG_ARRAY_pix(size,alpha_label,e,variables);
			else if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_NG_FnPix_pix(size,alpha_label,e,variables);*/
		}
		
		e -> minimize();
	
		//<-- Byung 5/14 Debugging
		for (int ii(0); ii < m_num_labels; ii++){
			m_classExistFlag[ii] = false;
		}
		m_classExistFlag[0] = true; // void always exists
		//--> Byung
	

		for ( i = 0; i < num; i++ )
		{
			if ( m_labeling[pixels[i]] != alpha_label )
			{
				if ( e->get_var(variables[i]) == 0 ){
					m_labeling[pixels[i]] = alpha_label;
					m_classExistFlag[alpha_label] = true;
				}
			}
			m_lookupPixVar[pixels[i]] = -1;
		}

		//<-- Byung 5/14 Debugging
		for ( i = 0; i < num; i++ )
		{
			m_classExistFlag[m_labeling[i]] = true;
		}
		//--> Byung

		delete [] variables;
	}

	delete e;

	return(compute_energy());
}

/**************************************************************************************/

void GCoptimization::add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label) //Msun: 3/5 add weights, add m_dataloss
{
	for ( int i = 0; i < size; i++ ){
		EnergyTermType total_cost_0 = (m_datacost(m_lookupPixVar[i],alpha_label)*m_seg_weight)+m_losscost(m_lookupPixVar[i],alpha_label);
		EnergyTermType total_cost_1 = (m_datacost(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]])*m_seg_weight)+ m_losscost(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]);
		// in case total_cost_0 or 1 are negative values.
		EnergyTermType bias = 0.0;
		if (total_cost_0<bias)
			bias = total_cost_0;
		if (total_cost_1<bias)
			bias = total_cost_1;
		
		e -> add_term1(variables[i], total_cost_0-bias,
		                             total_cost_1-bias);
	}
	
}

/**************************************************************************************/

void GCoptimization::add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label) //Msun: 3/5 add weights
{
	for ( int i = 0; i < size; i++ )
		e -> add_term1(variables[i], m_dataFnPix(m_lookupPixVar[i],alpha_label)*m_seg_weight,
		                             m_dataFnPix(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]])*m_seg_weight);

}
/**************************************************************************************/

void GCoptimization::add_t_links_FnCoord(Energy *e,Energy::Var *variables,int size,LabelType alpha_label) //Msun: 3/5 add weights
{
	int x,y;

	for ( int i = 0; i < size; i++ )
	{

		y = m_lookupPixVar[i]/m_width;
		x = m_lookupPixVar[i] - y*m_width;

		e -> add_term1(variables[i], m_dataFnCoord(x,y,alpha_label)*m_seg_weight,
		                             m_dataFnCoord(x,y,m_labeling[m_lookupPixVar[i]])*m_seg_weight);
	}

}

/**************************************************************************************/
// <!--Msun
/**************************************************************************************/

bool GCoptimization::CheckRegularity( int det_bool1, int det_bool2, int PtInd){

	//printf("Check regularity.\n"); fflush(stdout);
	bool regularFlag = false;
	if ( m_det_bool_negation[ det_bool1] == m_det_bool_negation[ det_bool2]){
		//printf("Cond1 entered.\n"); fflush(stdout);
		if( (m_PairThingWeights[ PtInd][0]+m_PairThingWeights[ PtInd][3])<=
		    (m_PairThingWeights[ PtInd][1]+m_PairThingWeights[ PtInd][2]))
			regularFlag = true;
	}else{
		//printf("Cond2 entered.\n"); fflush(stdout);
		if( (m_PairThingWeights[ PtInd][0]+m_PairThingWeights[ PtInd][3])>=
		    (m_PairThingWeights[ PtInd][1]+m_PairThingWeights[ PtInd][2]))
			regularFlag = true;
	}
	return regularFlag;
}

GCoptimization::PixelType GCoptimization::SetDetBoolNegLookUp( LabelType alpha_label){

	int ii,jj,PairThingInd,Pick;

	// simple set the negation of det_bool
	//printf("m_nHigherDet = %d\n", m_nHigherDet); fflush(stdout);
	for ( ii =0; ii < m_nHigherDet ; ii ++){
		if (!m_det_fix_bool[ii])
			m_det_bool_lookupValid[ii]=0; // default set to zero
		else
			m_det_bool_lookupValid[ii]=-1; // det_indicator should be fixed // Msun 9/13 added

		if (m_det_labeling[ii] == alpha_label)
			m_det_bool_negation[ii] =true;
		else
			m_det_bool_negation[ii] =false;
	}

	// decide which det_bool variables to fix. Now simply use all of them.
	//printf("Decide which det_bool among %d variables to fix. ",m_nPairThing); fflush(stdout);
	for ( ii=0; ii< m_nPairThing; ii++){
		//printf("ii=%d\n",ii); fflush(stdout);
		PairThingInd = m_PairThingTable[ii];
		//mexPrintf("PairThingInd=%d ",PairThingInd);
		//printf("PairThingInd=%d\n",PairThingInd); fflush(stdout);
		//printf("m_PairThingIndice[ PairThingInd ][0]=%d\n", m_PairThingIndice[ PairThingInd ][0]); fflush(stdout);
		//printf("m_PairThingIndice[ PairThingInd ][1]=%d\n", m_PairThingIndice[ PairThingInd ][1]); fflush(stdout);
		//printf("m_det_bool_lookupValid[ m_PairThingIndice[ PairThingInd ][0] ] = %d\n", m_det_bool_lookupValid[ m_PairThingIndice[ PairThingInd ][0] ]); fflush(stdout);
		if (m_det_bool_lookupValid[ m_PairThingIndice[ PairThingInd ][0] ]<0 ||
			m_det_bool_lookupValid[ m_PairThingIndice[ PairThingInd ][1] ]<0) // if any of the det_bool is fixed, then no need to fix the other
		{ 
			//printf("NO need to fix the other.\n"); fflush(stdout);
			continue;
		}
		
		//printf("Need to fix the other.\n"); fflush(stdout);


		if (!CheckRegularity(m_PairThingIndice[ PairThingInd][0], m_PairThingIndice[PairThingInd][1], PairThingInd)){ // if not regular, randomly fix one det_bool
			//printf("not regular.\n"); fflush(stdout);
#ifdef RAND_COMPILE
      		Pick = rand()%2;
#else
      		Pick = 0;
#endif 
			m_det_bool_lookupValid[ m_PairThingIndice[ PairThingInd ][Pick] ] = -1;
			//printf("m_PairThingIndice[%d][%d]=%d\n",PairThingInd,Pick,m_PairThingIndice[ PairThingInd ][Pick]); fflush(stdout);
		}
	}

	//printf("Rewrite m_det_bool_lookupValid...\n"); fflush(stdout);
	PixelType size = 0;
	for ( ii =0; ii < m_nHigherDet ; ii ++){
		if (m_det_bool_lookupValid[ii]>=0){
			m_det_bool_lookupValid[ii]=size;
			size++;
		}
	}
	return size;
}

void GCoptimization::perform_alpha_expansion_QPBO(LabelType alpha_label)
{	
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("in perform_alpha_expansion_QPBO()\n");
#endif
#else
	printf("in perform_alpha_expansion_QPBO()\n"); fflush(stdout);
#endif
	PixelType i,size = 0; 
	
	for ( i = 0; i < m_num_pixels; i++ )
	{
		if ( m_labeling[i] != alpha_label )
		{
			m_lookupPixVar[size] = i;
			m_lookupValidPix[i] = size;
			size++; // count number of valid pixel (i.e., the label does not equal to alpha_label)
		}else{
			m_lookupValidPix[i] = -1;
		}
	}
	
	int num_pairs;
	if ( m_grid_graph )
		num_pairs = count_neighbor_Grid( size, alpha_label);
	else
		num_pairs = count_neighbor_NGrid( size, alpha_label);
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("num_pairs=%d\n",num_pairs);
#endif
#else
	printf("num_pairs=%d\n",num_pairs); fflush(stdout);
#endif
	QPBO<REAL> *e = new QPBO<REAL>(size, num_pairs); // Msun: to-do get xxx: the number of edges
		
	if ( size > 0 ) 
	{
		e->AddNode(size);
		int *variables = (int *) new int[size+m_nHigher*2]; // Msun: each HOP add at most 2 more variables

		for ( i = 0; i < size; i++ )
			variables[i] = i; // initialize nodes in the graph

		/* add unary terms, see also m_datacost, m_dataFnPix, m_dataFnCoord---->*/
		if ( m_dataType == ARRAY ) add_t_links_ARRAY_QPBO(e,variables,size,alpha_label);
		else  if  ( m_dataType == FUNCTION_PIX ) add_t_links_FnPix_QPBO(e,variables,size,alpha_label);
		else  add_t_links_FnCoord_QPBO(e,variables,size,alpha_label);
		/* <----*/

		/* add higher-order-terms */
		if (m_hopType == 1){

			set_up_expansion_energy_QPBO_HOP( size, alpha_label, e, variables);
		}else{
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("Using DET formulation");
#endif
#else
			printf("Using DET formulation"); fflush(stdout);
#endif
			set_up_expansion_energy_QPBO_DET( size, alpha_label, e, variables);
		}
		/* <----*/
	
		/* add pair-wise smooth terms, see also m_smoothcost, m_smoothFnPix---->*/ //Msun: To-Do add higher order term
		if ( m_grid_graph )
		{
			if ( m_smoothType == ARRAY )
			{
				if (m_varying_weights) set_up_expansion_energy_G_ARRAY_VW_QPBO(size,alpha_label,e,variables);
				else set_up_expansion_energy_G_ARRAY_QPBO(size,alpha_label,e,variables);
			}
			else if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_G_FnPix_QPBO(size,alpha_label,e,variables);
			else  set_up_expansion_energy_G_FnCoord_QPBO(size,alpha_label,e,variables);
			
		}
		else
		{
			if ( m_smoothType == ARRAY ) set_up_expansion_energy_NG_ARRAY_QPBO(size,alpha_label,e,variables);
			else if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_NG_FnPix_QPBO(size,alpha_label,e,variables);
		}
		/* <----*/
	
		e -> MergeParallelEdges();// merge duplicated pair-wise terms
		e -> Solve();// solve graph cut
		e ->ComputeWeakPersistencies();
		REAL eng = e ->ComputeTwiceEnergy();
#ifdef MEX_COMPILE
#ifdef LOG_ON
		mexPrintf("QPBO eng=%f\n",eng/2.);
#endif
#else
		printf("QPBO eng=%f\n",eng/2.);
#endif
	
		// get result, modify m_labeling
		for ( i = 0,size = 0; i < m_num_pixels; i++ )
		{
			if ( m_labeling[i] != alpha_label )
			{
				if ( e->GetLabel(variables[size]) == 0 ){
					m_labeling[i] = alpha_label;
					m_classExistFlag[ alpha_label] = true;
				}

				size++;
			}
		}

		delete [] variables;
	}

	delete e;
}
/**************************************************************************************/
// for general HOP
void GCoptimization::set_up_expansion_energy_HOP(int size, LabelType alpha_label, Energy *e, Energy::Var *variables) //Msun: 3/5 added weights
{
	EnergyTermType lambda_a, lambda_d, lambda_m, gamma_d, number_old;
	LabelType maxLabel;
	int i, j;
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("m_nHigher=%d\n",m_nHigher);
#endif
#else
	printf("m_nHigher=%d\n",m_nHigher); fflush(stdout);
#endif
	for(i = 0;i < m_nHigher; i++)
	{
		maxLabel = getMaxLabel(i); // get dominant label  % Msun: to-do
		
		lambda_m = m_higherCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		lambda_a = m_higherCost[i * (m_num_labels + 1) + alpha_label];  // gamma_label (alpha)
		//if ((lambda_m - lambda_a) <0) return;

		variables[2*i+size] = e->add_variable(); // add auxilary node m_1 
		e->add_term1( variables[2 * i + size],0,(lambda_m - lambda_a)*m_seg_weight); // r_1
		for(j = 0; j < m_higherElements[i]; j++)
		{
			//if (m_lookupValidPix[m_higherIndex[i][j]]>= size)
			//	mexPrintf("m_lookupValidPix=%d ",m_lookupValidPix[m_higherIndex[i][j]]);
			if (m_lookupValidPix[m_higherIndex[i][j]]>=0) // not equal to alpha_label
				e->add_term2( variables[2*i + size], 
					variables[ m_lookupValidPix[m_higherIndex[i][j]]], 0,
					(m_higherWeights[i][j]*(lambda_m - lambda_a) / m_higherTruncation[i])*m_seg_weight, 0, 0); 
		}   

		//mexPrintf("maxLabel=%d lambda_m=%f lambda_a=%f ",maxLabel,lambda_m,lambda_a);
		if((maxLabel != -1) && (maxLabel != alpha_label)) // there exist a dominant label
		{   
			number_old = cardinality(i, maxLabel);// weights of nodes labeld dominant in current potential (w_i influencing)
			gamma_d = m_higherCost[i * (m_num_labels + 1) + maxLabel];
			lambda_d = gamma_d + 
				(m_higherP[i] - number_old) // R_d including weights
				*( m_higherCost[i * (m_num_labels + 1) + m_num_labels] - m_higherCost[i * (m_num_labels + 1) + maxLabel]) // gamma_max - gamma_d
				*(1 / m_higherTruncation[i]); // 1/Q


			if ((lambda_m - lambda_d) <0) return;
			variables[2*i+size+1] = e->add_variable(); // auxilary node m_0
			e->add_term1(variables[2 * i + size + 1],(lambda_m - lambda_d)*m_seg_weight,0); //weight r_0
			for(j = 0; j < m_higherElements[i]; j++)
				if (m_labeling[m_higherIndex[i][j]] == maxLabel) // connect dominant-labeled nodes to m_0
					e->add_term2(variables[2*i + size+1 ], variables[ m_lookupValidPix[m_higherIndex[i][j]] ], 
						0 , 0, (m_higherWeights[i][j]*(lambda_m - gamma_d) / m_higherTruncation[i])*m_seg_weight, 0); 
			
			//mexPrintf("number_old=%f gamma_d=%f lambda_d=%f \n", number_old,gamma_d,lambda_d);
		}
	}
}

// for detection HOP
void GCoptimization::set_up_expansion_energy_Det_HOP2(LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables) //Msun: 3/5 add weights
{
	EnergyTermType gamma_m, neg_det_bool_w, thing_uw;
	int i, j;
//	EnergyTermType tmp_Energy; // Byung Added for debug
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("m_nHigherDet=%d\n",m_nHigherDet);
#endif
#else
	printf("m_nHigherDet=%d\n",m_nHigherDet); fflush(stdout);
#endif
	for(i = 0;i < m_nHigherDet; i++)
	{
/*// <-- Byung Added		
//GCoptimization::EnergyType GCoptimization::giveDetHOP2Energy() // Byung Added for debug
		tmp_Energy = giveDetHOP2Energy();
		if (tmp_Energy != 0) {
			printf("Break:Debug\n");
		}

// Byung Added	-->*/


		gamma_m = m_higherDetCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		thing_uw = m_higherDetUW[i];
		EnergyTermType tmpE = m_higherDetP[i]-m_higherDetTruncation[i];
		if (!m_det_bool_negation[i]){ // alpha \not= det_laveling & negation is false
			if (m_det_bool_lookupValid[i]<0){ // det_bool is fixed
				if ( m_det_bool_label[i] == false){
					for(j = 0; j < m_higherDetElements[i]; j++){
						// penalty for y is false
						if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i] ){ // equal to det_label
							e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
								(m_higherDetWeights[i][j]*gamma_m/tmpE)*m_det_class_weights[m_det_labeling[i]]); // add element unary potential
						}
					}
				}else{
					for(j = 0; j < m_higherDetElements[i]; j++){
						// penalty for y is true
						if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i] ){ // equal to det_label
							e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 
								(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0); // add element unary potential
						}
					}
				}
			}else{ // det_bool not fixed
				neg_det_bool_w = m_higherDetP[i];
				for(j = 0; j < m_higherDetElements[i]; j++)
				{
					if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i] ){ // equal to det_label
						neg_det_bool_w -= m_higherDetWeights[i][j];
						// penalty for y is true
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 0,
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0); 
						// penalty for y is false
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
							(m_higherDetWeights[i][j]*gamma_m/tmpE)*m_det_class_weights[m_det_labeling[i]], 0, 0); 
					}
				}
#ifdef NO_UNARY_DET_COMPILE
				e->add_term1( det_variables[ m_det_bool_lookupValid[i] ], gamma_m*m_det_class_weights[m_det_labeling[i]],
					(neg_det_bool_w*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]]);
#else
				e->add_term1( det_variables[ m_det_bool_lookupValid[i] ], thing_uw*m_det_class_weights[m_det_labeling[i]],
					(neg_det_bool_w*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]]);
#endif
			}
		}else{ // alpha == det_labeling
			if (m_det_bool_lookupValid[i]<0){
				if ( m_det_bool_label[i] == false){
					// penalty for y is false
					for(j = 0; j < m_higherDetElements[i]; j++){
						if (m_labeling[ m_higherDetIndex[i][j]] != m_det_labeling[i] ){ // equal to det_label
							e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]],
								( m_higherDetWeights[i][j]*gamma_m/tmpE)*m_det_class_weights[m_det_labeling[i]],0 ); // add element unary potential
						}
					}
				}else{
					// penalty for y is true
					for(j = 0; j < m_higherDetElements[i]; j++){
						if (m_labeling[ m_higherDetIndex[i][j]] != m_det_labeling[i] ){ // equal to det_label
							e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
								(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]]); // add element unary potential
						}
					}
				}
			}else{
				neg_det_bool_w = m_higherDetP[i];
				for(j = 0; j < m_higherDetElements[i]; j++)
				{
					if (m_labeling[ m_higherDetIndex[i][j]] != m_det_labeling[i] ){ // not equal to det_label
						neg_det_bool_w -= m_higherDetWeights[i][j];
						// penalty for y is true
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0, 0); 
						// penalty for y is false
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 0, 
							(m_higherDetWeights[i][j]*gamma_m/tmpE)*m_det_class_weights[m_det_labeling[i]], 0); 
					}
				}   
#ifdef NO_UNARY_DET_COMPILE
				e->add_term1( det_variables[m_det_bool_lookupValid[i]], 0, (gamma_m+(neg_det_bool_w*gamma_m/tmpE))*m_det_class_weights[m_det_labeling[i]] );
#else
				e->add_term1( det_variables[m_det_bool_lookupValid[i]], 0, (thing_uw+(neg_det_bool_w*gamma_m/tmpE))*m_det_class_weights[m_det_labeling[i]] );
#endif
			}
		}
	}
}

// for detection HOP
void GCoptimization::set_up_expansion_energy_Det_HOP(LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables) //Msun: 3/5 added wiehgts
{
	EnergyTermType gamma_m, neg_det_bool_w;
	int i, j;

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("m_nHigherDet=%d\n",m_nHigherDet);
#endif
#else
	printf("m_nHigherDet=%d\n",m_nHigherDet); fflush(stdout);
#endif
	for(i = 0;i < m_nHigherDet; i++)
	{
		gamma_m = m_higherDetCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		if (!m_det_bool_negation[i]){ // alpha \not= det_laveling & negation is false
			if (m_det_bool_lookupValid[i]<0){ // det_bool is fixed
				if ( m_det_bool_label[i] == false)
					continue;
				for(j = 0; j < m_higherDetElements[i]; j++){
					if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i] ){ // equal to det_label
						e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0); // add element unary potential
					}
				}
			}else{ // det_bool not fixed
				neg_det_bool_w = m_higherDetP[i];
				for(j = 0; j < m_higherDetElements[i]; j++)
				{
					if (m_labeling[ m_higherDetIndex[i][j]] == m_det_labeling[i] ){ // equal to det_label
						neg_det_bool_w -= m_higherDetWeights[i][j];
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 0,
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0); 
					}
				}   
				e->add_term1( det_variables[ m_det_bool_lookupValid[i] ], gamma_m*m_det_class_weights[m_det_labeling[i]],
					(neg_det_bool_w*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]]);
			}
		}else{ // alpha == det_labeling
			if (m_det_bool_lookupValid[i]<0){
				if ( m_det_bool_label[i] == false)
					continue;
				for(j = 0; j < m_higherDetElements[i]; j++){
					if (m_labeling[ m_higherDetIndex[i][j]] != m_det_labeling[i] ){ // equal to det_label
						e->add_term1( variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]]); // add element unary potential
					}
				}
			}else{
				for(j = 0; j < m_higherDetElements[i]; j++)
				{
					if (m_labeling[ m_higherDetIndex[i][j]] != m_det_labeling[i] ){ // not equal to det_label
						e->add_term2( det_variables[ m_det_bool_lookupValid[i] ], 
							variables[ m_lookupValidPix[m_higherDetIndex[i][j]]], 0, 
							(m_higherDetWeights[i][j]*gamma_m/m_higherDetTruncation[i])*m_det_class_weights[m_det_labeling[i]], 0, 0); 
					}
				}   
				e->add_term1( det_variables[m_det_bool_lookupValid[i]], 0, gamma_m*m_det_class_weights[m_det_labeling[i]]);
			}
		}
	}
}

// for detection HOP
void GCoptimization::set_up_expansion_energy_Occ_HOP(LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables)
{
	EnergyTermType gamma_m, neg_det_bool_w;
	int i, j;

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("m_nHigherOcc=%d\n",m_nHigherOcc);
#endif
#else
	printf("m_nHigherOcc=%d\n",m_nHigherOcc); fflush(stdout);
#endif
	for(i = 0;i < m_nHigherOcc; i++)
	{
		gamma_m = m_higherOccCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		if (!m_det_bool_negation[ m_Occl2HigherDetIndex[i]]){ // alpha \not= det_laveling
			if (m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ]<0){ // det_bool is fixed
				if ( m_det_bool_label[ m_Occl2HigherDetIndex[i]] == false)
					continue;
				for(j = 0; j < m_higherOccElements[i]; j++){
					if (m_labeling[ m_higherOccIndex[i][j]] == m_det_labeling[ m_Occl2HigherDetIndex[i]] ){ // equal to det_label
						e->add_term1( variables[ m_lookupValidPix[m_higherOccIndex[i][j]]], m_higherOccWeights[i][j]*gamma_m/m_higherOccTruncation[i], 0); // add element unary potential
					}
				}
			}else{ // det_bool not fixed
				neg_det_bool_w = m_higherOccP[i];
				for(j = 0; j < m_higherOccElements[i]; j++)
				{
					if (m_labeling[ m_higherOccIndex[i][j]] == m_det_labeling[ m_Occl2HigherDetIndex[i] ] ){ // equal to det_label
						neg_det_bool_w -= m_higherOccWeights[i][j];
						e->add_term2( det_variables[ m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ] ], 
							variables[ m_lookupValidPix[m_higherOccIndex[i][j]]], 0, 0,
							m_higherOccWeights[i][j]*gamma_m/m_higherOccTruncation[i], 0); 
					}
				}   
				e->add_term1( det_variables[ m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ] ], gamma_m,
					(neg_det_bool_w*gamma_m/m_higherOccTruncation[i]));
			}
		}else{
			if (m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ]<0){
				if ( m_det_bool_label[ m_Occl2HigherDetIndex[i]] == false)
					continue;
				for(j = 0; j < m_higherOccElements[i]; j++){
					if (m_labeling[ m_higherOccIndex[i][j]] != m_det_labeling[m_Occl2HigherDetIndex[i]] ){ // equal to det_label
						e->add_term1( variables[ m_lookupValidPix[m_higherOccIndex[i][j]]], 0, m_higherOccWeights[i][j]*gamma_m/m_higherOccTruncation[i]); // add element unary potential
					}
				}
			}else{
				for(j = 0; j < m_higherOccElements[i]; j++)
				{
					if (m_labeling[ m_higherOccIndex[i][j]] != m_det_labeling[m_Occl2HigherDetIndex[i]] ){ // not equal to det_label
						e->add_term2( det_variables[ m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ] ], 
							variables[ m_lookupValidPix[m_higherOccIndex[i][j]]], 0, 
							m_higherOccWeights[i][j]*gamma_m/m_higherOccTruncation[i], 0, 0); 
					}
				}   
				e->add_term1( det_variables[m_det_bool_lookupValid[ m_Occl2HigherDetIndex[i] ]], 0, gamma_m);
			}
		}
	}
}

// pair-wise thing indicator potential
void GCoptimization::set_up_expansion_energy_Pair_Thing(LabelType alpha_label, Energy *e, Energy::Var *det_variables, int det_size)
{
	int ii, pt1, pt2;
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("m_nPairThing=%d\n",m_nPairThing);
#endif
#else
	printf("m_nPairThing=%d\n",m_nPairThing); fflush(stdout);
#endif
	for (ii=0; ii <m_nPairThing; ii++)
	{
		pt1 = m_PairThingIndice[ii][0];
		pt2 = m_PairThingIndice[ii][1];
		//if (pt1 >= m_nHigherDet || pt2 >= m_nHigherDet){
		//	mexPrintf("pt1=%d,pt2=%d,m_nHigherDet=%d\n",pt1,pt2,m_nHigherDet);
		//	terminateOnError(true, "in set_up_expansion_energy_Pair_Thing pt1 pt2 error");
		//}
		//if (m_det_bool_lookupValid[pt1]>= det_size || m_det_bool_lookupValid[pt2]>= det_size){
		//	mexPrintf("det_bool_lkup1=%d,det_bool_lkup2=%d,det_size=%d\n",m_det_bool_lookupValid[pt1],m_det_bool_lookupValid[pt2],det_size);
		//	terminateOnError(true, "in set_up_expansion_energy_Pair_Thing m_det_bool_lookupValid error");
		//}
		// find min negative bias ====== 2/27 new code ======
		GCoptimization::EnergyType bias = 0;
//		for (int kk=0; kk<4; kk++){
//			if (m_PairThingWeights[ii][kk]< bias){
//				bias = m_PairThingWeights[ii][kk];
//			}
//		}
		// END find min negative bias ====== 2/27 new code ==

		if ( m_det_bool_lookupValid[pt1]>=0 && m_det_bool_lookupValid[pt2]>=0){ // should set pair-wise term
			if ( !m_det_bool_negation[pt1] && !m_det_bool_negation[pt2]){
				e->add_term2( det_variables[m_det_bool_lookupValid[pt1]],
						      det_variables[m_det_bool_lookupValid[pt2]], 
							  m_PairThingWeights[ii][0]-bias, m_PairThingWeights[ii][1]-bias,
							  m_PairThingWeights[ii][2]-bias, m_PairThingWeights[ii][3]-bias);
			}else if ( m_det_bool_negation[pt1] && !m_det_bool_negation[pt2]){
				e->add_term2( det_variables[m_det_bool_lookupValid[pt1]],
						      det_variables[m_det_bool_lookupValid[pt2]], 
							  m_PairThingWeights[ii][2]-bias, m_PairThingWeights[ii][3]-bias,
							  m_PairThingWeights[ii][0]-bias, m_PairThingWeights[ii][1]-bias);
			}else if ( !m_det_bool_negation[pt1] && m_det_bool_negation[pt2]){
				e->add_term2( det_variables[m_det_bool_lookupValid[pt1]],
						      det_variables[m_det_bool_lookupValid[pt2]], 
							  m_PairThingWeights[ii][1]-bias, m_PairThingWeights[ii][0]-bias,
							  m_PairThingWeights[ii][3]-bias, m_PairThingWeights[ii][2]-bias);
			}else{
				e->add_term2( det_variables[m_det_bool_lookupValid[pt1]],
						      det_variables[m_det_bool_lookupValid[pt2]], 
							  m_PairThingWeights[ii][3]-bias, m_PairThingWeights[ii][2]-bias,
							  m_PairThingWeights[ii][1]-bias, m_PairThingWeights[ii][0]-bias);
			}
		}else if ( m_det_bool_lookupValid[pt1]<0 && m_det_bool_lookupValid[pt2]>=0){ // pt1 is fixed
			if ( !m_det_bool_negation[pt2]){
				if ( m_det_bool_label[pt1]){
					e->add_term1( det_variables[m_det_bool_lookupValid[pt2]],
						m_PairThingWeights[ii][2]-bias, m_PairThingWeights[ii][3]-bias);
				}else{
					e->add_term1( det_variables[m_det_bool_lookupValid[pt2]],
						m_PairThingWeights[ii][0]-bias, m_PairThingWeights[ii][1]-bias);
				}
			}else{
				if ( m_det_bool_label[pt1]){
					e->add_term1( det_variables[m_det_bool_lookupValid[pt2]],
						m_PairThingWeights[ii][3]-bias, m_PairThingWeights[ii][2]-bias);
				}else{
					e->add_term1( det_variables[m_det_bool_lookupValid[pt2]],
						m_PairThingWeights[ii][1]-bias, m_PairThingWeights[ii][0]-bias);
				}
			}
		}else if ( m_det_bool_lookupValid[pt1]>=0 && m_det_bool_lookupValid[pt2]<0){ // pt2 is fixed
			if ( !m_det_bool_negation[pt1]){
				if ( m_det_bool_label[pt2]){
					e->add_term1( det_variables[m_det_bool_lookupValid[pt1]],
						m_PairThingWeights[ii][1]-bias, m_PairThingWeights[ii][3]-bias);
				}else{
					e->add_term1( det_variables[m_det_bool_lookupValid[pt1]],
						m_PairThingWeights[ii][0]-bias, m_PairThingWeights[ii][2]-bias);
				}
			}else{
				if ( m_det_bool_label[pt2]){
					e->add_term1( det_variables[m_det_bool_lookupValid[pt1]],
						m_PairThingWeights[ii][3]-bias, m_PairThingWeights[ii][1]-bias);
				}else{
					e->add_term1( det_variables[m_det_bool_lookupValid[pt1]],
						m_PairThingWeights[ii][2]-bias, m_PairThingWeights[ii][0]-bias);
				}
			}
		}
	}
}

// for Class Co
void GCoptimization::set_up_expansion_energy_Class_Co(LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *co_variables)
{
	bool switchFlag = true;

	if (!m_useClassCoFlag)
		return;

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("in GCoptimization::set_up_expansion_energy_Class_Co");
#endif
#endif
	int ii,jj;
	EnergyTermType kappa_alpha  = m_class_co_u_w[ alpha_label];
	EnergyTermType kappa_l[m_num_labels];
	for (ii=0; ii< m_num_labels; ii++) kappa_l[ii] = 0;

	// loop through all lables
	for (ii=0; ii< m_num_labels; ii++){
		if (ii != alpha_label && m_classExistFlag[ii]){ //get kappa_l and set z_l unart
			kappa_alpha += m_class_co_p_w[alpha_label+ii*m_num_labels];
			for (jj=0; jj< m_num_labels; jj++){
				if (jj!= alpha_label && ii!=jj && m_classExistFlag[jj]){
					if (kappa_l[ii]> m_class_co_p_w[ii+jj*m_num_labels])
						kappa_l[ii] = m_class_co_p_w[ii+jj*m_num_labels];
				}
			}
		}
	}

	// get unary term
	if (switchFlag){ // t=0, means expand to alpha
		e->add_term1( co_variables[alpha_label], kappa_alpha, 0);
	}else{
		e->add_term1( co_variables[alpha_label], 0, kappa_alpha);
	}
	for (ii=0; ii< m_num_labels; ii++){
		if (ii != alpha_label && m_classExistFlag[ii]){ //get kappa_l and set z_l unart
			if (switchFlag){ // t=0, means expand to alpha
				e->add_term1( co_variables[ii], 0, kappa_l[ii]);
			}else{
				e->add_term1( co_variables[ii], kappa_l[ii], 0);
			}
		}
	}

//	if (false){	
	// loop through all elements
	for (ii=0; ii< m_num_pixels; ii++){
		if (m_lookupValidPix[ii]>=0){// set pair-wise term
			if (switchFlag){ // t=0, means expand to alpha
				e->add_term2( co_variables[alpha_label], variables[m_lookupValidPix[ii]],
					0, 0, kappa_alpha, 0);
				e->add_term2( co_variables[m_labeling[ii]], variables[m_lookupValidPix[ii]],
					0, kappa_l[m_labeling[ii]], 0, 0);
			}else{ //t=1, means expand to alpha
				e->add_term2( co_variables[alpha_label], variables[m_lookupValidPix[ii]],
					0, kappa_alpha, 0, 0);
				e->add_term2( co_variables[m_labeling[ii]], variables[m_lookupValidPix[ii]],
					0, 0, kappa_l[m_labeling[ii]], 0);
			}
		}
	}
//	}
}

// for Class Co Compressed
void GCoptimization::set_up_expansion_energy_CompClass_Co(LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *co_variables)
{
	bool switchFlag = true;

	if (!m_useClassCoFlag || m_num_comp_labels==0)
		return;

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("in GCoptimization::set_up_expansion_energy_CompClass_Co");
#endif
#endif
	int ii,jj,sub_ii,sub_jj;
	LabelType comp_alpha_label = m_label2comp_label[ alpha_label]; // get the compressed label
	EnergyTermType kappa_alpha  = m_class_co_u_w[ comp_alpha_label];
	EnergyTermType kappa_l[m_num_comp_labels];
	for (ii=0; ii< m_num_comp_labels; ii++) kappa_l[ii] = 0;

	// loop through all lables
	for (ii=0; ii< m_num_comp_labels; ii++){
		for (sub_ii=0; sub_ii< m_num_dup_labels[ii]; sub_ii++){
			LabelType target_label_ii = m_comp_label2label[ii][sub_ii];
			if (ii != comp_alpha_label && m_classExistFlag[ target_label_ii ] ){ //get kappa_l and set z_l unart
				kappa_alpha += m_class_co_p_w[comp_alpha_label+ii*m_num_comp_labels];
				for (jj=0; jj< m_num_comp_labels; jj++){
					for (sub_jj=0; sub_jj< m_num_dup_labels[jj]; sub_jj++){
						LabelType target_label_jj = m_comp_label2label[jj][sub_jj];
						if (jj!= comp_alpha_label && ii!=jj && m_classExistFlag[target_label_jj] ){
							if (kappa_l[ii]> m_class_co_p_w[ii+jj*m_num_comp_labels])
								kappa_l[ii] = m_class_co_p_w[ii+jj*m_num_comp_labels];
							
							break;
						}
					}
				}
				break;
			}
		}
	}

	// get unary term
	if (switchFlag){ // t=0, means expand to alpha
		e->add_term1( co_variables[comp_alpha_label], kappa_alpha, 0);
	}else{
		e->add_term1( co_variables[comp_alpha_label], 0, kappa_alpha);
	}
	for (ii=0; ii< m_num_comp_labels; ii++){
		for (sub_ii=0; sub_ii< m_num_dup_labels[ii]; sub_ii++){
			LabelType target_label_ii = m_comp_label2label[ii][sub_ii];
			if (ii != comp_alpha_label && m_classExistFlag[target_label_ii]){ //get kappa_l and set z_l unart
				if (switchFlag){ // t=0, means expand to alpha
					e->add_term1( co_variables[ii], 0, kappa_l[ii]);
				}else{
					e->add_term1( co_variables[ii], kappa_l[ii], 0);
				}
				break;
			}
		}
	}

	// get pair-wise term
	// loop through all elements
	for (ii=0; ii< m_num_pixels; ii++){
		if (m_lookupValidPix[ii]>=0){// set pair-wise term
			LabelType comp_m_labeling = m_label2comp_label[ m_labeling[ii]]; // get the compressed label
			if (switchFlag){ // t=0, means expand to alpha
				e->add_term2( co_variables[comp_alpha_label], variables[m_lookupValidPix[ii]],
					0, 0, kappa_alpha, 0);
				e->add_term2( co_variables[comp_m_labeling], variables[m_lookupValidPix[ii]],
					0, kappa_l[comp_m_labeling], 0, 0);
			}else{ //t=1, means expand to alpha
				e->add_term2( co_variables[comp_alpha_label], variables[m_lookupValidPix[ii]],
					0, kappa_alpha, 0, 0);
				e->add_term2( co_variables[comp_m_labeling], variables[m_lookupValidPix[ii]],
					0, 0, kappa_l[comp_m_labeling], 0);
			}
		}
	}
}

//--------------------- QPBO ----------------------
void GCoptimization::set_up_expansion_energy_QPBO_HOP(int size, LabelType alpha_label, QPBO<REAL> *e, int *variables)
{
	EnergyTermType lambda_a, lambda_d, lambda_m, gamma_d, number_old;
	LabelType maxLabel;
	int i, j, count;
	count = 0;

	//mexPrintf("m_nHigher=%d\n",m_nHigher);
	for(i = 0;i < m_nHigher; i++)
	{
		maxLabel = getMaxLabel(i); // get dominant label  % Msun: to-do
	
		lambda_m = m_higherCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		lambda_a = m_higherCost[i * (m_num_labels + 1) + alpha_label];  // gamma_label (alpha)

		e->AddNode(1); // add auxilary node m_1 
		variables[count+size] = count+size;
		// combine add_tweights add_edge to AddPairwiseTerm
		e->AddUnaryTerm( variables[count + size],0,lambda_m - lambda_a); // r_1
		for(j = 0; j < m_higherElements[i]; j++)
		{   
			if (m_lookupValidPix[m_higherIndex[i][j]]>=0) // not equal to alpha_label
			{
				e->AddPairwiseTerm( variables[count + size], 
					variables[ m_lookupValidPix[m_higherIndex[i][j]]], 0,
					m_higherWeights[i][j]*(lambda_m - lambda_a) / m_higherTruncation[i], 0, 0); 
			}
		}   
		count++;

		if((maxLabel != -1) && (maxLabel != alpha_label)) // there exist a dominant label
		{   
			number_old = cardinality(i, maxLabel);// weights of nodes labeld dominant in current potential (w_i influencing)
			gamma_d = m_higherCost[i * (m_num_labels + 1) + maxLabel];
			lambda_d = gamma_d + 
				(m_higherP[i] - number_old) // R_d including weights
				*( m_higherCost[i * (m_num_labels + 1) + m_num_labels] - m_higherCost[i * (m_num_labels + 1) + maxLabel]) // gamma_max - gamma_d
				*(1 / m_higherTruncation[i]); // 1/Q


			e->AddNode(1); // auxilary node m_0
			variables[count+size] = count+size;
			e->AddUnaryTerm(variables[count + size ],lambda_m - lambda_d,0); //weight r_0
			for(j = 0; j < m_higherElements[i]; j++)
				if (m_labeling[m_higherIndex[i][j]] == maxLabel) // connect dominant-labeled nodes to m_0
					e->AddPairwiseTerm(variables[count + size ], variables[ m_lookupValidPix[m_higherIndex[i][j]] ], 
						0 , 0, m_higherWeights[i][j]*(lambda_m - gamma_d) / m_higherTruncation[i], 0); 
			count++;
		}
	}
	//mexPrintf("2*m_nHigher=%d count=%d\n",2*m_nHigher,count);
}

void GCoptimization::set_up_expansion_energy_QPBO_DET(int size, LabelType alpha_label, QPBO<REAL> *e, int *variables)
{
	EnergyTermType lambda_a, lambda_d, lambda_m, gamma_d, number_old;
	LabelType maxLabel;
	int i, j, count;
	count = 0;

	//mexPrintf("m_nHigher=%d\n",m_nHigher);
	for(i = 0;i < m_nHigher; i++)
	{
		maxLabel = getMaxLabel(i); // get dominant label  % Msun: to-do
	
		lambda_m = m_higherCost[i * (m_num_labels + 1) + m_num_labels]; // gamma_max 
		lambda_a = m_higherCost[i * (m_num_labels + 1) + alpha_label];  // gamma_label (alpha)

		e->AddNode(1); // add auxilary node m_0, which is the indicator for detection
		variables[count+size] = count+size;
		// combine add_tweights add_edge to AddPairwiseTerm
		e->AddUnaryTerm( variables[count + size],lambda_m - lambda_a,0); // r_1
		for(j = 0; j < m_higherElements[i]; j++)
		{   
			if (m_lookupValidPix[m_higherIndex[i][j]]>=0) // not equal to alpha_label
				//e->AddPairwiseTerm( variables[count + size], 
				//	variables[ m_lookupValidPix[m_higherIndex[i][j]]], 0,
				//	m_higherWeights[i][j]*(lambda_m - lambda_a) / m_higherTruncation[i], 0, 0); 
				e->AddPairwiseTerm( variables[count + size], 
					variables[ m_lookupValidPix[m_higherIndex[i][j]]], 0,0,
					m_higherWeights[i][j]*(lambda_m - lambda_a) / m_higherTruncation[i], 0); 
		}   

		if((maxLabel != -1) && (maxLabel != alpha_label)) // there exist a dominant label
		{   
			number_old = cardinality(i, maxLabel);// weights of nodes labeld dominant in current potential (w_i influencing)
			gamma_d = m_higherCost[i * (m_num_labels + 1) + maxLabel];
			lambda_d = gamma_d + 
				(m_higherP[i] - number_old) // R_d including weights
				*( m_higherCost[i * (m_num_labels + 1) + m_num_labels] - m_higherCost[i * (m_num_labels + 1) + maxLabel]) // gamma_max - gamma_d
				*(1 / m_higherTruncation[i]); // 1/Q


			e->AddUnaryTerm(variables[count + size ],lambda_m - lambda_d,0); //weight r_0
			for(j = 0; j < m_higherElements[i]; j++)
				if (m_labeling[m_higherIndex[i][j]] == maxLabel) // connect dominant-labeled nodes to m_0
					e->AddPairwiseTerm(variables[count + size ], variables[ m_lookupValidPix[m_higherIndex[i][j]] ], 
						0 , 0, m_higherWeights[i][j]*(lambda_m - gamma_d) / m_higherTruncation[i], 0); 
		}
		count++;

	}
	//mexPrintf("2*m_nHigher=%d count=%d\n",2*m_nHigher,count);
}

/*
* For HOpotential i, choose dominant label (-1 if there is not dominant label) - can be at most one
* label d s.t.: W(c_d) > P - Q_d,  
*/
int GCoptimization::getMaxLabel(int i)
{
	int j;
	EnergyTermType *num_labels = new EnergyTermType[m_num_labels];

	for(j = 0;j < m_num_labels; j++)
		num_labels[j] = 0;

	for(j = 0;j < m_higherElements[i]; j++)
		num_labels[m_labeling[m_higherIndex[i][j]]]+= m_higherWeights[i][j];

	EnergyTermType number = 0;
	int maxLabel;

	for(j = 0;j < m_num_labels; j++)
	{
		if(number <= num_labels[j])
		{
			number = num_labels[j];
			maxLabel = j;
		}
	}

	delete[] num_labels;
	if(number > (m_higherP[i] - m_higherTruncation[i])) // Assumes same Q for all labels
		return maxLabel; 
	else 
	return -1;
}
/*
* For HO-potential i sum w_j delta_label(x_j)
*/
GCoptimization::EnergyTermType GCoptimization::cardinality(int i, int label)
{
	int j;
	EnergyTermType count_label = 0;

	for(j = 0;j<m_higherElements[i]; j++)
		if(m_labeling[m_higherIndex[i][j]] == label)  
			count_label+=m_higherWeights[i][j];
	
	return count_label;
}

/**************************************************************************************/

void GCoptimization::add_t_links_ARRAY_QPBO(QPBO<REAL> *e, int *variables,int size,LabelType alpha_label)
{
	for ( int i = 0; i < size; i++ )
		e -> AddUnaryTerm( variables[i], m_datacost(m_lookupPixVar[i],alpha_label),
		                             m_datacost(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}

/**************************************************************************************/

void GCoptimization::add_t_links_FnPix_QPBO(QPBO<REAL> *e, int *variables,int size,LabelType alpha_label)
{
	for ( int i = 0; i < size; i++ )
		e -> AddUnaryTerm(variables[i], m_dataFnPix(m_lookupPixVar[i],alpha_label),
		                             m_dataFnPix(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}
/**************************************************************************************/

void GCoptimization::add_t_links_FnCoord_QPBO(QPBO<REAL> *e,int *variables,int size,LabelType alpha_label)
{
	int x,y;

	for ( int i = 0; i < size; i++ )
	{

		y = m_lookupPixVar[i]/m_width;
		x = m_lookupPixVar[i] - y*m_width;

		e -> AddUnaryTerm(variables[i], m_dataFnCoord(x,y,alpha_label),
		                             m_dataFnCoord(x,y,m_labeling[m_lookupPixVar[i]]));
	}

}

/**************************************************************************************/

/**********************************************************************************************/
/* Performs alpha-expansion for non regular grid graph for case when energy terms are NOT     */
/* specified by a function */

void GCoptimization::set_up_expansion_energy_NG_ARRAY_QPBO(int size, LabelType alpha_label,QPBO<REAL> *e,int *variables )
{
	EnergyTermType weight;
	Neighbor *tmp;
	int i,nPix,pix;;

	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		m_lookupPixVar[pix] = i;

		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront();
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				weight = tmp->weight;
				
				if ( m_labeling[nPix] != alpha_label )
				{
					if ( pix < nPix )
						e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothcost(alpha_label,alpha_label)*weight,
									  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
									  m_smoothcost(m_labeling[pix],alpha_label)*weight,
									  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
				}
				else
					e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,alpha_label)*weight,
  					                       m_smoothcost(m_labeling[pix],alpha_label)*weight);
				
			}
		}
	}

	
}

/**********************************************************************************************/
/* Performs alpha-expansion for non regular grid graph for case when energy terms are        */
/* specified by a function */

void GCoptimization::set_up_expansion_energy_NG_FnPix_QPBO(int size, LabelType alpha_label, QPBO<REAL> *e,int *variables )
{
	Neighbor *tmp;
	int i,nPix,pix;
	


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		m_lookupPixVar[pix] = i;


		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront(); //??
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				
				if ( m_labeling[nPix] != alpha_label )
				{
					if ( pix < nPix )
						e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
									  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
									  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label),
									  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]));
				}
				else
					e ->AddUnaryTerm(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
						                       m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
				
			}
		}
	}
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY_VW_QPBO(int size, LabelType alpha_label,QPBO<REAL> *e,
													  int *variables )
{
	int i,nPix,pix,x,y, EdgeId;
	EnergyTermType weight;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
			weight = m_horizWeights[pix];
			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
			else   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
			weight = m_vertWeights[pix];
			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight ,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight ,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight );
			else   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[nPix],
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_horizWeights[nPix]);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,alpha_label)*m_vertWeights[nPix],
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_vertWeights[nPix]);
		}	
			
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("EdgeId=%d\n",EdgeId);	
#endif
#else
	printf("EdgeId=%d\n",EdgeId);	
#endif
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY_QPBO(int size, LabelType alpha_label, QPBO<REAL> *e,
													 int *variables )
{
	int i,nPix,pix,x,y,EdgeId;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label),
							  m_smoothcost(alpha_label,m_labeling[nPix]),
							  m_smoothcost(m_labeling[pix],alpha_label),
							  m_smoothcost(m_labeling[pix],m_labeling[nPix]));
			else   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
				                 m_smoothcost(m_labeling[pix],alpha_label));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label) ,
							  m_smoothcost(alpha_label,m_labeling[nPix]),
							  m_smoothcost(m_labeling[pix],alpha_label) ,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
			else   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
				                 m_smoothcost(m_labeling[pix],alpha_label));
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
			   	                 m_smoothcost(m_labeling[pix],alpha_label) );
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothcost(alpha_label,alpha_label),
			   	                 m_smoothcost(m_labeling[pix],alpha_label));
		}	
			
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("EdgeId=%d\n",EdgeId);	
#endif
#else
	printf("EdgeId=%d\n",EdgeId);	
#endif
	
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_FnPix_QPBO(int size, LabelType alpha_label, QPBO<REAL> *e,
													 int *variables )
{
	int i,nPix,pix,x,y, EdgeId;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
							  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
							  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label),
							  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]));
			else   e ->AddUnaryTerm(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
				                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothFnPix(pix,nPix,alpha_label,alpha_label) ,
							  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
							  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label) ,
							  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]) );
			else   e ->AddUnaryTerm(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
				                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
			   	                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label) );
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
			   	                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
		}	
			
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("EdgeId=%d\n",EdgeId);	
#endif
#else
	printf("EdgeId=%d\n",EdgeId);	
#endif

}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_FnCoord_QPBO(int size, LabelType alpha_label,QPBO<REAL> *e,
													 int *variables )
{
	int i,nPix,pix,x,y, EdgeId;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_horz_cost(x,y,alpha_label,alpha_label),
							  m_horz_cost(x,y,alpha_label,m_labeling[nPix]),
							  m_horz_cost(x,y,m_labeling[pix],alpha_label),
							  m_horz_cost(x,y,m_labeling[pix],m_labeling[nPix]));
			else   e ->AddUnaryTerm(variables[i],m_horz_cost(x,y,alpha_label,alpha_label),
				                 m_horz_cost(x,y,m_labeling[pix],alpha_label));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				EdgeId = e ->AddPairwiseTerm(variables[i],variables[m_lookupPixVar[nPix]],
					          m_vert_cost(x,y,alpha_label,alpha_label) ,
							  m_vert_cost(x,y,alpha_label,m_labeling[nPix]),
							  m_vert_cost(x,y,m_labeling[pix],alpha_label) ,
							  m_vert_cost(x,y,m_labeling[pix],m_labeling[nPix]) );
			else   e ->AddUnaryTerm(variables[i],m_vert_cost(x,y,alpha_label,alpha_label),
				                 m_vert_cost(x,y,m_labeling[pix],alpha_label));
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_horz_cost(x-1,y,alpha_label,m_labeling[nPix]),
			   	                 m_horz_cost(x-1,y,m_labeling[pix],alpha_label) );
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->AddUnaryTerm(variables[i],m_vert_cost(x,y-1,alpha_label,alpha_label),
			   	                 m_vert_cost(x,y-1,m_labeling[pix],alpha_label));
		}	
	}
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("EdgeId=%d\n",EdgeId);	
#endif
#else
	printf("EdgeId=%d\n",EdgeId);	
#endif

}

/**************************************************************************************/
int GCoptimization::count_neighbor_Grid(int size, LabelType alpha_label)
{
	int i,nPix,pix,x,y,count;
	count = 0;
	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
			if ( m_labeling[nPix] != alpha_label )
				count ++;
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
			if ( m_labeling[nPix] != alpha_label )
				count ++;
		}	
	}

	// count hop edges
	for(i = 0;i < m_nHigher; i++)	
	{
		count = count+2*m_higherElements[i];
	}   
	
	return count;
}
/**************************************************************************************/
int GCoptimization::count_neighbor_NGrid(int size, LabelType alpha_label)
{
	Neighbor *tmp;
	int i,nPix,pix,count;
	count = 0;

	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		//m_lookupPixVar[pix] = i;

		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront();
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				
				if ( m_labeling[nPix] != alpha_label )
				{
					if ( pix < nPix )
						count ++;
				}
			}
		}
	}

	// count hop edges
	for(i = 0;i < m_nHigher; i++)	
	{
		count = count+2*m_higherElements[i];
	}
   
	return count;
}
/**************************************************************************************/
// Msun -->

/**************************************************************************************/
void GCoptimization::perform_alpha_expansion(LabelType alpha_label)
{
	PixelType i, det_size, size = 0;
#ifdef MEX_COMPILE
	Energy *e = new Energy(mexErrMsgTxt);
#else
	Energy *e = new Energy();
#endif
#ifdef LOG_ON
	printf("Starting func perform_alpha_expansion...\n"); fflush(stdout);
#endif
	for ( i = 0; i < m_num_pixels; i++ )
	{
		//printf("(%d/%d)", i, m_num_pixels);
		if ( m_labeling[i] != alpha_label )
		{
			//printf("m_labeling[%d]=%d != alpha_label=%d\n", i, m_labeling[i],  alpha_label);
			m_lookupPixVar[size] = i;
			m_lookupValidPix[i] = size;
			size++; // count number of valid pixel (i.e., the label does not equal to alpha_label)
		}else{
			//printf("m_labeling[%d]=%d == alpha_label=%d\n", i, m_labeling[i],  alpha_label);
			m_lookupValidPix[i] = -1;
		}
	}
	
#ifdef LOG_ON
	printf("Starting SetDetBoolNegLookup...\n"); fflush(stdout);
#endif
	// <--Msun: add setting m_det_bool_negation and m_det_bool_lookupValid
	det_size = SetDetBoolNegLookUp(alpha_label);
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("det_size=%d\n",det_size);
#endif
#else
	printf("det_size=%d\n",det_size); fflush(stdout);
#endif
	// --> Msun
	
	if ( size > 0 ) 
	{
		Energy::Var *variables = (Energy::Var *) new Energy::Var[size+2*m_nHigher];

		for ( i = 0; i < size; i++ )
			variables[i] = e ->add_variable(); // initialize nodes in the graph

		/* Msun initialize det_variable */
		Energy::Var *det_variables = (Energy::Var *) new Energy::Var[det_size];
		for ( i = 0; i < det_size; i++ )
			det_variables[i] = e ->add_variable(); // initialize nodes in the graph
	
		/* Msun initilize ClassCo variable */
		Energy::Var *co_variables = (Energy::Var *) new Energy::Var[m_num_labels];// at most m_num_labels
		if (m_num_comp_labels==0){
			for ( i = 0; i < m_num_labels; i++ )
				co_variables[i] = e ->add_variable(); // initialize nodes in the graph
		}else{
			for ( i = 0; i < m_num_comp_labels; i++ )
				co_variables[i] = e ->add_variable(); // initialize nodes in the graph
		}
	
		/* add unary terms, see also m_datacost, m_dataFnPix, m_dataFnCoord---->*/
		if ( m_dataType == ARRAY ) add_t_links_ARRAY(e,variables,size,alpha_label);
		else  if  ( m_dataType == FUNCTION_PIX ) add_t_links_FnPix(e,variables,size,alpha_label);
		else  add_t_links_FnCoord(e,variables,size,alpha_label);
		/* <----*/

		/* Msun: add Hop terms */
		set_up_expansion_energy_HOP( size, alpha_label, e, variables);
#ifdef LOG_ON
		printf("Energy setting(HOP) is done.\n"); fflush(stdout);
#endif
		/* Msun: add Det Hop terms */
#ifdef OLD_DET_COMPILE
		set_up_expansion_energy_Det_HOP( alpha_label, e, variables, det_variables);
#else
		set_up_expansion_energy_Det_HOP2( alpha_label, e, variables, det_variables);
#endif
#ifdef LOG_ON
		printf("Energy setting(Det_HOP) is done.\n"); fflush(stdout);
#endif
		/* Msun: add Occ Hop terms */
		set_up_expansion_energy_Occ_HOP( alpha_label, e, variables, det_variables);
#ifdef LOG_ON
		printf("Energy setting(Occ_HOP) is done.\n"); fflush(stdout);
#endif
		/* Msun: add Pair Thing terms */
		set_up_expansion_energy_Pair_Thing( alpha_label, e, det_variables, det_size);
#ifdef LOG_ON
		printf("Energy setting(Pair_Thing) is done.\n"); fflush(stdout);
#endif
		/* Msun: add Class Co term */
		if (m_num_comp_labels==0){
			set_up_expansion_energy_Class_Co( alpha_label, e, variables, co_variables);
		}else{
			set_up_expansion_energy_CompClass_Co( alpha_label, e, variables, co_variables);
		}


		/* add pair-wise smooth terms, see also m_smoothcost, m_smoothFnPix---->*/ //Msun: To-Do add higher order term
		if ( m_grid_graph )
		{
			if ( m_smoothType == ARRAY )
			{
				if (m_varying_weights) set_up_expansion_energy_G_ARRAY_VW(size,alpha_label,e,variables);
				else set_up_expansion_energy_G_ARRAY(size,alpha_label,e,variables);
			}
			else if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_G_FnPix(size,alpha_label,e,variables);
			else  set_up_expansion_energy_G_FnCoord(size,alpha_label,e,variables);
			
		}
		else
		{
			if ( m_smoothType == ARRAY ) set_up_expansion_energy_NG_ARRAY(size,alpha_label,e,variables);
			else if ( m_smoothType == FUNCTION_PIX ) set_up_expansion_energy_NG_FnPix(size,alpha_label,e,variables);
		}
		/* <----*/
		
#ifdef LOG_ON
		printf("Energy setting is done.\n"); fflush(stdout);
#endif
		Energy::TotalValue Emin = e -> minimize();// solve graph cut

#ifdef LOG_ON
		printf("Graphcut is solved.\n"); fflush(stdout);
#endif

#ifdef MEX_COMPILE
		//mexPrintf("GCeng=%f\n",Emin);
#else
		//printf("GCeng=%f\n",Emin);
#endif


		//<-- Byung 5/14 Debugging
		for (int ii(0); ii < m_num_labels; ii++){
			m_classExistFlag[ii] = false;
		}
		m_classExistFlag[0] = true; // void always exists


		//--> Byung
		// get result, modify m_labeling
		for ( i = 0,size = 0; i < m_num_pixels; i++ )
		{
			if ( m_labeling[i] != alpha_label )
			{
				if ( e->get_var(variables[size]) == 0 ){
					m_labeling[i] = alpha_label;
					m_classExistFlag[alpha_label] = true;
				}
				size++;
			}
		}

		//<-- Byung 5/14 Debugging
		for ( i = 0; i < m_num_pixels; i++ )
		{
			m_classExistFlag[m_labeling[i]] = true;
		}
		//--> Byung



#ifdef LOG_ON
		printf("m_labeling is modified.\n"); fflush(stdout);
#endif
		// Msun: 2/2 added get indicator results
		for (i = 0; i< m_nHigherDet; i++){
			if ( m_det_bool_lookupValid[i]>=0){
				if (m_det_bool_negation[i]){
					if (e->get_var( det_variables[m_det_bool_lookupValid[i]])==1){
						m_det_bool_label[i] = false;
					}else{
						m_det_bool_label[i] = true;
					}
				}else{
					if (e->get_var( det_variables[m_det_bool_lookupValid[i]])==1){
						m_det_bool_label[i] = true;
					}else{
						m_det_bool_label[i] = false;
					}
				}
			}
		}

#ifdef LOG_ON
		printf("m_det_bool_label is modified.\n"); fflush(stdout);
#endif
		delete [] det_variables;
		delete [] co_variables;
		delete [] variables;
		//if (det_size !=0) delete [] det_variables;
		//delete [] co_variables;
		
	}

	delete e;
}

/**********************************************************************************************/
/* Performs alpha-expansion for non regular grid graph for case when energy terms are NOT     */
/* specified by a function */

void GCoptimization::set_up_expansion_energy_NG_ARRAY(int size, LabelType alpha_label,Energy *e,Energy::Var *variables ) //Msun: 3/5 added weights
{
	EnergyTermType weight;
	Neighbor *tmp;
	int i,nPix,pix;;



	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		m_lookupPixVar[pix] = i;

		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront();
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				weight = tmp->weight*m_seg_weight;
				
				if ( m_labeling[nPix] != alpha_label )
				{
					if ( pix < nPix )
						e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothcost(alpha_label,alpha_label)*weight,
									  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
									  m_smoothcost(m_labeling[pix],alpha_label)*weight,
									  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
				}
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*weight,
  					                       m_smoothcost(m_labeling[pix],alpha_label)*weight);
				
			}
		}
	}

	
}

/**********************************************************************************************/
/* Performs alpha-expansion for non regular grid graph for case when energy terms are        */
/* specified by a function */

void GCoptimization::set_up_expansion_energy_NG_FnPix(int size, LabelType alpha_label,Energy *e,Energy::Var *variables ) //Msun: 3/5 ad weights
{
	Neighbor *tmp;
	int i,nPix,pix;
	


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		m_lookupPixVar[pix] = i;


		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront(); //??
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				
				if ( m_labeling[nPix] != alpha_label )
				{
					if ( pix < nPix )
						e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothFnPix(pix,nPix,alpha_label,alpha_label)*m_seg_weight,
									  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
									  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight,
									  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix])*m_seg_weight);
				}
				else
					e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
						                       m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight);
				
			}
		}
	}
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY_VW(int size, LabelType alpha_label,Energy *e, //Msun: 3/5 ad weights
													  Energy::Var *variables )
{
	int i,nPix,pix,x,y;
	EnergyTermType weight;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
			weight = m_horizWeights[pix]*m_seg_weight;
			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
			weight = m_vertWeights[pix]*m_seg_weight;
			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight ,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight ,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight );
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[nPix]*m_seg_weight,
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_horizWeights[nPix]*m_seg_weight);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*m_vertWeights[nPix]*m_seg_weight,
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_vertWeights[nPix]*m_seg_weight);
		}	
			
	}
	
}





/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy *e, // Msun: 3/5 add weights
													 Energy::Var *variables )
{
	int i,nPix,pix,x,y;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*m_seg_weight,
							  m_smoothcost(alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_seg_weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*m_seg_weight,
							  m_smoothcost(alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_seg_weight,
				                 m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_seg_weight,
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*m_seg_weight,
			   	                 m_smoothcost(m_labeling[pix],alpha_label)*m_seg_weight);
		}	
			
	}
	
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY_VW_pix(int size, LabelType alpha_label,Energy *e,
													  Energy::Var *variables, PixelType *pixels, int num )
{
	int i,nPix,pix,x,y;
	EnergyTermType weight;


	for ( i = 0; i < num; i++ )
	{
		pix = pixels[i];
		if ( m_labeling[pix]!= alpha_label )
		{
			y = pix/m_width;
			x = pix - y*m_width;

			if ( x < m_width - 1 )
			{
				nPix = pix + 1;
				weight = m_horizWeights[pix];
				if ( m_labeling[nPix] != alpha_label && m_lookupPixVar[nPix] != -1 )
					e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
				else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
			}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
			weight = m_vertWeights[pix];
			if ( m_labeling[nPix] != alpha_label && m_lookupPixVar[nPix] != -1 )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label)*weight ,
							  m_smoothcost(alpha_label,m_labeling[nPix])*weight,
							  m_smoothcost(m_labeling[pix],alpha_label)*weight ,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight );
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
				                 m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label || m_lookupPixVar[nPix] == -1)
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[nPix],
			   	                 m_smoothcost(m_labeling[pix],m_labeling[nPix])*m_horizWeights[nPix]);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label ||	m_lookupPixVar[nPix] == -1)
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*m_vertWeights[nPix],
			   	                 m_smoothcost(m_labeling[pix],m_labeling[nPix])*m_vertWeights[nPix]);
		}	
		}	
	}
	
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_ARRAY_pix(int size, LabelType alpha_label,Energy *e,
													 Energy::Var *variables, PixelType *pixels,
													 int num)
{
	int i,nPix,pix,x,y;


	for ( i = 0; i < num; i++ )
	{
		pix = pixels[i];
		y = pix/m_width;
		x = pix - y*m_width;


		if ( m_labeling[pix]!= alpha_label )
		{
			if ( x < m_width - 1 )
			{
				nPix = pix + 1;

				if ( m_labeling[nPix] != alpha_label && m_lookupPixVar[pix] != -1)
					e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label),
							  m_smoothcost(alpha_label,m_labeling[nPix]),
							  m_smoothcost(m_labeling[pix],alpha_label),
							  m_smoothcost(m_labeling[pix],m_labeling[nPix]));
				else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
				                 m_smoothcost(m_labeling[pix],m_labeling[nPix]));
			}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label && m_lookupPixVar[pix] != -1)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothcost(alpha_label,alpha_label) ,
							  m_smoothcost(alpha_label,m_labeling[nPix]),
							  m_smoothcost(m_labeling[pix],alpha_label) ,
							  m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
			else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
				                 m_smoothcost(m_labeling[pix],m_labeling[nPix]));
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label || m_lookupPixVar[nPix] == -1)
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
			   	                 m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label || m_lookupPixVar[nPix] == -1)
			   e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),
			   	                 m_smoothcost(m_labeling[pix],m_labeling[nPix]));
		}	
			
	}
	}
}


/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_FnPix(int size, LabelType alpha_label,Energy *e, // Msun: 3/5 add weights
													 Energy::Var *variables )
{
	int i,nPix,pix,x,y;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothFnPix(pix,nPix,alpha_label,alpha_label)*m_seg_weight,
							  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight,
							  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
				                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_smoothFnPix(pix,nPix,alpha_label,alpha_label)*m_seg_weight,
							  m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight,
							  m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
				                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix])*m_seg_weight,
			   	                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,alpha_label)*m_seg_weight,
			   	                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label)*m_seg_weight);
		}	
			
	}
	
}

/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void GCoptimization::set_up_expansion_energy_G_FnCoord(int size, LabelType alpha_label,Energy *e, //Msun:3/5 added weights
													 Energy::Var *variables )
{
	int i,nPix,pix,x,y;


	for ( i = size - 1; i >= 0; i-- )
	{
		pix = m_lookupPixVar[i];
		y = pix/m_width;
		x = pix - y*m_width;

		m_lookupPixVar[pix] = i;

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_horz_cost(x,y,alpha_label,alpha_label)*m_seg_weight,
							  m_horz_cost(x,y,alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_horz_cost(x,y,m_labeling[pix],alpha_label)*m_seg_weight,
							  m_horz_cost(x,y,m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_horz_cost(x,y,alpha_label,alpha_label)*m_seg_weight,
				                 m_horz_cost(x,y,m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;

			if ( m_labeling[nPix] != alpha_label )
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
					          m_vert_cost(x,y,alpha_label,alpha_label)*m_seg_weight,
							  m_vert_cost(x,y,alpha_label,m_labeling[nPix])*m_seg_weight,
							  m_vert_cost(x,y,m_labeling[pix],alpha_label)*m_seg_weight,
							  m_vert_cost(x,y,m_labeling[pix],m_labeling[nPix])*m_seg_weight);
			else   e ->add_term1(variables[i],m_vert_cost(x,y,alpha_label,alpha_label)*m_seg_weight,
				                 m_vert_cost(x,y,m_labeling[pix],alpha_label)*m_seg_weight);
		}	
		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_horz_cost(x-1,y,alpha_label,m_labeling[nPix])*m_seg_weight,
			   	                 m_horz_cost(x-1,y,m_labeling[pix],alpha_label)*m_seg_weight);
		}	

		if ( y > 0 )
		{
			nPix = pix - m_width;
	
			if ( m_labeling[nPix] == alpha_label )
			   e ->add_term1(variables[i],m_vert_cost(x,y-1,alpha_label,alpha_label)*m_seg_weight,
			   	                 m_vert_cost(x,y-1,m_labeling[pix],alpha_label)*m_seg_weight);
		}	
	}
}



/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::oneExpansionIteration()
{
	int next;

	terminateOnError( m_dataType == NONE,"You have to set up the data cost before running optimization");
	terminateOnError( m_smoothType == NONE,"You have to set up the smoothness cost before running optimization");


	if (m_random_label_order) scramble_label_table();

	if (m_random_pair_thing_order) scramble_pair_thing_table();
	

#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("in oneExpansionIteration");
#endif
#else
	printf("in oneExpansionIteration (m_num_labels=%d)", m_num_labels); fflush(stdout);
#endif
	for (next = 0;  next < m_num_labels;  next++ ){
		EnergyType new_energy,old_energy;
	
		//printf("Before calc energ\n"); fflush(stdout);
		old_energy = compute_energy();

		//perform_alpha_expansion(m_labelTable[next]);
		//printf("Before expansion starts\n"); fflush(stdout);
		if (m_solver == GCsolver) perform_alpha_expansion(m_labelTable[next]);
		else if (m_solver == QPBOsolver) {
		perform_alpha_expansion_QPBO(m_labelTable[next]);
		}
		else terminateOnError(1,"Did not specify the solver!");

		new_energy = compute_energy();
		if (old_energy<new_energy){
#ifdef MEX_COMPILE
#ifdef LOG_ON
			mexPrintf("Label=%d,EngWWrapOld=%f,New=%f\n",m_labelTable[next],old_energy,new_energy);
#endif
#else
			printf("Label(%d)=%d,EngWWrapOld=%f,New=%f\n",next,m_labelTable[next],old_energy,new_energy);
#endif
		}			
    }
#ifdef MEX_COMPILE
#ifdef LOG_ON
	mexPrintf("out oneExpansionIteration");
#endif
#else
	printf("out oneExpansionIteration");
#endif
	

	
	return(compute_energy());
}

/**************************************************************************************/

void GCoptimization::setNeighbors(PixelType pixel1, int pixel2, EnergyTermType weight)
{

	assert(pixel1 < m_num_pixels && pixel1 >= 0 && pixel2 < m_num_pixels && pixel2 >= 0);
	assert(m_grid_graph == 0);

	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	temp1->weight  = weight;
	temp1->to_node = pixel2;

	temp2->weight  = weight;
	temp2->to_node = pixel1;

	m_neighbors[pixel1].addFront(temp1);
	m_neighbors[pixel2].addFront(temp2);

	// <-- Byung
	//delete temp1;
	//delete temp2;
	// Byung -->
	
}

/**************************************************************************************/

void GCoptimization::setNeighbors(PixelType pixel1, int pixel2)
{

	assert(pixel1 < m_num_pixels && pixel1 >= 0 && pixel2 < m_num_pixels && pixel2 >= 0);
	assert(m_grid_graph == 0);
	

	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	temp1->weight  = (EnergyTermType) 1;
	temp1->to_node = pixel2;

	temp2->weight  = (EnergyTermType) 1;
	temp2->to_node = pixel1;

	m_neighbors[pixel1].addFront(temp1);
	m_neighbors[pixel2].addFront(temp2);

	// <-- Byung
	//delete temp1;
	//delete temp2;
	// Byung -->
	
}

/**************************************************************************************/

GCoptimization::~GCoptimization()
{

    // <!-- bagon
    // mexWarnMsgTxt("Calling destructor"); /* bagon added */
    if ( m_dataType == ARRAY ){
        delete [] m_datacost;
        delete [] m_losscost;
	}
    if ( m_smoothType == ARRAY )
        delete [] m_smoothcost;
    if ( m_varying_weights == 1 ) {
        delete [] m_vertWeights;
        delete [] m_horizWeights;
    }
    // bagon -->
	
   
	// <--Muns: 
	int ii;
	if (m_nHigher !=0){
		delete [] m_lookupValidPix;
		delete [] m_higherCost;
		delete [] m_higherTruncation;
		delete [] m_higherP;
		delete [] m_higherElements;
		for (ii=0;ii<m_nHigher;ii++){
			delete [] m_higherWeights[ii];
			delete [] m_higherIndex[ii];
		}
		delete [] m_higherWeights;
		delete [] m_higherIndex;
	}


	if (m_nHigherDet !=0){
		delete [] m_det_bool_label;
		delete [] m_det_bool_negation;
		delete [] m_det_bool_lookupValid;
		delete [] m_higherDetCost;
		delete [] m_higherDetUW;
		delete [] m_det_labeling;
		delete [] m_higherDetTruncation;
		delete [] m_higherDetP;
		delete [] m_higherDetElements;
		for (ii=0;ii<m_nHigherDet;ii++){
			delete [] m_higherDetWeights[ii];
			delete [] m_higherDetIndex[ii];
		}
		delete [] m_higherDetWeights;
		delete [] m_higherDetIndex;
		delete [] m_det_class_weights;
	}	


	if (m_nHigherOcc !=0){
		delete [] m_Occl2HigherDetIndex;
		delete [] m_higherOccCost;
		delete [] m_higherOccTruncation;
		delete [] m_higherOccP;
		delete [] m_higherOccElements;
		for (ii=0;ii<m_nHigherOcc;ii++){
			delete [] m_higherOccWeights[ii];
			delete [] m_higherOccIndex[ii];
		}
		delete [] m_higherOccWeights;
		delete [] m_higherOccIndex;
	}


	if (m_nPairThing !=0){
		delete [] m_PairThingTable;
		for (ii=0;ii<m_nPairThing;ii++){
			delete [] m_PairThingIndice[ii];
			delete [] m_PairThingWeights[ii];
		}
		delete [] m_PairThingIndice;
		delete [] m_PairThingWeights;
	}


	if (m_useClassCoFlag){
		delete [] m_class_co_u_w;
		delete [] m_class_co_p_w;
		delete [] m_classExistFlag;
	}

	if (m_num_comp_labels !=0){
		for (ii=0;ii<m_num_comp_labels;ii++){
			delete [] m_comp_label2label[ii];
		}
		delete [] m_comp_label2label;
		delete [] m_label2comp_label;
		delete [] m_num_dup_labels;
	}
	// Msun -->
 
	if ( m_deleteLabeling ) 
		delete [] m_labeling;


	if ( m_dataInput == SET_INDIVIDUALLY ){
		delete [] m_datacost;
		delete [] m_losscost;
	}

	if ( m_smoothInput == SET_INDIVIDUALLY )
			delete [] m_smoothcost;
 		
	if ( ! m_grid_graph )
		delete [] m_neighbors;			


	delete [] m_labelTable;
	delete [] m_lookupPixVar;

			
}


/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::swap(int max_num_iterations)
{
	return(start_swap(max_num_iterations)); 
}

/**************************************************************************************/


GCoptimization::EnergyType GCoptimization::swap()
{
	return(start_swap(MAX_INTT));
}

/**************************************************************************************/


GCoptimization::EnergyType GCoptimization::start_swap(int max_num_iterations )
{
	
	int curr_cycle = 1;
	EnergyType new_energy,old_energy;
	

	new_energy = compute_energy();

	old_energy = (new_energy+1)*10; // BAGON changed init value to exceed current energy by factor of 10 (thanks to A. Khan)    

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
		new_energy = oneSwapIteration();
		
		curr_cycle++;	
	}

	return(new_energy);
}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::oneSwapIteration()
{
	
	int next,next1;
   
	if (m_random_label_order) scramble_label_table();
		

	for (next = 0;  next < m_num_labels;  next++ )
		for (next1 = m_num_labels - 1;  next1 >= 0;  next1-- )
			if ( m_labelTable[next] < m_labelTable[next1] )
				perform_alpha_beta_swap(m_labelTable[next],m_labelTable[next1]);

	return(compute_energy());
}

/**************************************************************************************/

GCoptimization::EnergyType GCoptimization::alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
	terminateOnError( alpha_label < 0 || alpha_label >= m_num_labels || beta_label < 0 || beta_label >= m_num_labels,
		"Illegal Label to Expand On");
	perform_alpha_beta_swap(alpha_label,beta_label);
	return(compute_energy());
}
/**************************************************************************************/

void GCoptimization::add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size,
											LabelType alpha_label, LabelType beta_label,
											PixelType *pixels)
{
	for ( int i = 0; i < size; i++ )
		e -> add_term1(variables[i], m_datacost(pixels[i],alpha_label),
									 m_datacost(pixels[i],beta_label));

}
	
/**************************************************************************************/

void GCoptimization::add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size,
											LabelType alpha_label, LabelType beta_label,
											PixelType *pixels)
{
	for ( int i = 0; i < size; i++ )
		e -> add_term1(variables[i], m_dataFnPix(pixels[i],alpha_label),
									 m_dataFnPix(pixels[i],beta_label));

}
/**************************************************************************************/

void GCoptimization::add_t_links_FnCoord_swap(Energy *e,Energy::Var *variables,int size,
											  LabelType alpha_label, LabelType beta_label,
											  PixelType *pixels)
{
	int x,y;

	for ( int i = 0; i < size; i++ )
	{

		y = pixels[i]/m_width;
		x = pixels[i] - y*m_width;

		e -> add_term1(variables[i], m_dataFnCoord(x,y,alpha_label),m_dataFnCoord(x,y,beta_label));
	}
}


/**************************************************************************************/

void  GCoptimization::perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
	PixelType i,size = 0;
#ifdef MEX_COMPILE
	Energy *e = new Energy(mexErrMsgTxt);
#else
	Energy *e = new Energy();
#endif
	PixelType *pixels = new PixelType[m_num_pixels];
	

	for ( i = 0; i < m_num_pixels; i++ )
	{
		if ( m_labeling[i] == alpha_label || m_labeling[i] == beta_label)
		{
			pixels[size]    = i;
			m_lookupPixVar[i] = size;
			size++;
		}
	}

	if ( size == 0 )
	{
		delete e;
		delete [] pixels;
		return;
	}


	Energy::Var *variables = (Energy::Var *) new Energy::Var[size];


	for ( i = 0; i < size; i++ )
		variables[i] = e ->add_variable();

	if ( m_dataType == ARRAY ) add_t_links_ARRAY_swap(e,variables,size,alpha_label,beta_label,pixels);
	else  if  ( m_dataType == FUNCTION_PIX ) add_t_links_FnPix_swap(e,variables,size,alpha_label,beta_label,pixels);
	else  add_t_links_FnCoord_swap(e,variables,size,alpha_label,beta_label,pixels);



	if ( m_grid_graph )
	{
		if ( m_smoothType == ARRAY )
		{
			if (m_varying_weights) set_up_swap_energy_G_ARRAY_VW(size,alpha_label,beta_label,pixels,e,variables);
			else set_up_swap_energy_G_ARRAY(size,alpha_label,beta_label,pixels,e,variables);
		}
		else if ( m_smoothType == FUNCTION_PIX ) set_up_swap_energy_G_FnPix(size,alpha_label,beta_label,pixels,e,variables);
		else  set_up_swap_energy_G_FnCoord(size,alpha_label,beta_label,pixels,e,variables);
			
	}
	else
	{
		if ( m_smoothType == ARRAY ) set_up_swap_energy_NG_ARRAY(size,alpha_label,beta_label,pixels,e,variables);
		else if ( m_smoothType == FUNCTION_PIX ) set_up_swap_energy_NG_FnPix(size,alpha_label,beta_label,pixels,e,variables);
	}
		

	e -> minimize();

	//<-- Byung 5/14 Debugging
	for (int ii(0); ii < m_num_labels; ii++){
		m_classExistFlag[ii] = false;
	}
	m_classExistFlag[0] = true; // void always exists
	//--> Byung


	for ( i = 0; i < size; i++ )
		if ( e->get_var(variables[i]) == 0 ){
			m_labeling[pixels[i]] = alpha_label;
			m_classExistFlag[alpha_label] = true;
		}else{
			 m_labeling[pixels[i]] = beta_label;
			m_classExistFlag[beta_label] = true;
		}
	//<-- Byung 5/14 Debugging
	for ( i = 0; i < m_num_pixels; i++ )
	{
		m_classExistFlag[m_labeling[i]] = true;
	}
	//--> Byung


	delete [] variables;
	delete [] pixels;
	delete e;

}

/**************************************************************************************/

void GCoptimization::set_up_swap_energy_NG_ARRAY(int size,LabelType alpha_label,LabelType beta_label,
												 PixelType *pixels,Energy* e, Energy::Var *variables)
{
	PixelType nPix,pix,i;
	EnergyTermType weight;
	Neighbor *tmp;
	


	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront();
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				weight = tmp->weight;
				
				if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				{
					if ( pix < nPix )
						e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothcost(alpha_label,alpha_label)*weight,
									  m_smoothcost(alpha_label,beta_label)*weight,
									  m_smoothcost(beta_label,alpha_label)*weight,
									  m_smoothcost(beta_label,beta_label)*weight);
				}
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
						                       m_smoothcost(beta_label,m_labeling[nPix])*weight);
			}
		}
	}
}

/**************************************************************************************/

void GCoptimization::set_up_swap_energy_NG_FnPix(int size,LabelType alpha_label,LabelType beta_label,
												 PixelType *pixels,Energy* e, Energy::Var *variables)
{
	PixelType nPix,pix,i;
	Neighbor *tmp;
	

	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		if ( !m_neighbors[pix].isEmpty() )
		{
			m_neighbors[pix].setCursorFront();
			
			while ( m_neighbors[pix].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[pix].next());
				nPix   = tmp->to_node;
				
				if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				{
					if ( pix < nPix )
						e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
							          m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
									  m_smoothFnPix(pix,nPix,alpha_label,beta_label),
									  m_smoothFnPix(pix,nPix,beta_label,alpha_label),
									  m_smoothFnPix(pix,nPix,beta_label,beta_label) );
				}
				else
					e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
						                       m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
			}
		}
	}
}

/**************************************************************************************/

void GCoptimization::set_up_swap_energy_G_FnPix(int size,LabelType alpha_label,LabelType beta_label,
												PixelType *pixels,Energy* e, Energy::Var *variables)
{
	PixelType nPix,pix,i,x,y;

	
	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		y = pix/m_width;
		x = pix - y*m_width;

		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
							  m_smoothFnPix(pix,nPix,alpha_label,beta_label),
							  m_smoothFnPix(pix,nPix,beta_label,alpha_label),
							  m_smoothFnPix(pix,nPix,beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
 					                       m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
	
		}	
		if ( y > 0 )
		{
			nPix = pix - m_width;
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
							  m_smoothFnPix(pix,nPix,alpha_label,beta_label),
							  m_smoothFnPix(pix,nPix,beta_label,alpha_label),
							  m_smoothFnPix(pix,nPix,beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
 					                       m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
		}	

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
					e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
 					                           m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
				e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
				                           m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));

		}
	}
}
/**************************************************************************************/

void GCoptimization::set_up_swap_energy_G_FnCoord(int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels,
     											 Energy* e, Energy::Var *variables)
{  
	PixelType nPix,pix,i,x,y;

	
	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		y = pix/m_width;
		x = pix - y*m_width;

		if ( x > 0 )
		{
			nPix = pix - 1;
	
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_horz_cost(x-1,y,alpha_label,alpha_label),
							  m_horz_cost(x-1,y,alpha_label,beta_label),
							  m_horz_cost(x-1,y,beta_label,alpha_label),
							  m_horz_cost(x-1,y,beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_horz_cost(x-1,y,alpha_label,m_labeling[nPix]),
 					                           m_horz_cost(x-1,y,beta_label,m_labeling[nPix]));
	
		}	
		if ( y > 0 )
		{
			nPix = pix - m_width;
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_vert_cost(x,y-1,alpha_label,alpha_label),
							  m_vert_cost(x,y-1,alpha_label,beta_label),
							  m_vert_cost(x,y-1,beta_label,alpha_label),
							  m_vert_cost(x,y-1,beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_vert_cost(x,y-1,alpha_label,m_labeling[nPix]),
 					                       m_vert_cost(x,y-1,beta_label,m_labeling[nPix]));
		}	

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
					e ->add_term1(variables[i],m_horz_cost(x,y,alpha_label,m_labeling[nPix]),
 					                           m_horz_cost(x,y,beta_label,m_labeling[nPix]));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
				e ->add_term1(variables[i],m_vert_cost(x,y,alpha_label,m_labeling[nPix]),
				                           m_vert_cost(x,y,beta_label,m_labeling[nPix]));

		}
	}
}

/**************************************************************************************/

void GCoptimization::set_up_swap_energy_G_ARRAY_VW(int size,LabelType alpha_label,LabelType beta_label,
												   PixelType *pixels,Energy* e, Energy::Var *variables)
{
	PixelType nPix,pix,i,x,y;
	EnergyTermType weight;	



	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		y = pix/m_width;
		x = pix - y*m_width;

		if ( x > 0 )
		{
			nPix = pix - 1;
			weight = m_horizWeights[nPix];
	
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothcost(alpha_label,alpha_label)*weight,
							  m_smoothcost(alpha_label,beta_label)*weight,
							  m_smoothcost(beta_label,alpha_label)*weight,
							  m_smoothcost(beta_label,beta_label)*weight );
	
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
 					                       m_smoothcost(beta_label,m_labeling[nPix])*weight);
	
		}	
		if ( y > 0 )
		{
			nPix = pix - m_width;
			weight = m_vertWeights[nPix];

			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothcost(alpha_label,alpha_label)*weight,
							  m_smoothcost(alpha_label,beta_label)*weight,
							  m_smoothcost(beta_label,alpha_label)*weight,
							  m_smoothcost(beta_label,beta_label)*weight );
	
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
 					                       m_smoothcost(beta_label,m_labeling[nPix])*weight);
		}	

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[pix],
 					                           m_smoothcost(beta_label,m_labeling[nPix])*m_horizWeights[pix]);
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
				e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_vertWeights[pix],
				                           m_smoothcost(beta_label,m_labeling[nPix])*m_vertWeights[pix]);

		}
	}
}

/**************************************************************************************/

void GCoptimization::set_up_swap_energy_G_ARRAY(int size,LabelType alpha_label,LabelType beta_label,
											   PixelType *pixels,Energy* e, Energy::Var *variables)

{
	PixelType nPix,pix,i,x,y;
	


	for ( i = 0; i < size; i++ )
	{
		pix = pixels[i];
		y = pix/m_width;
		x = pix - y*m_width;

		if ( x > 0 )
		{
			nPix = pix - 1;
	
	
			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothcost(alpha_label,alpha_label),
							  m_smoothcost(alpha_label,beta_label),
							  m_smoothcost(beta_label,alpha_label),
							  m_smoothcost(beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
 					                       m_smoothcost(beta_label,m_labeling[nPix]));
	
		}	
		if ( y > 0 )
		{
			nPix = pix - m_width;

			if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
				e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
				              m_smoothcost(alpha_label,alpha_label),
							  m_smoothcost(alpha_label,beta_label),
							  m_smoothcost(beta_label,alpha_label),
							  m_smoothcost(beta_label,beta_label) );
	
				else
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
 					                       m_smoothcost(beta_label,m_labeling[nPix]));
		}	

		if ( x < m_width - 1 )
		{
			nPix = pix + 1;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
					e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
 					                           m_smoothcost(beta_label,m_labeling[nPix]));
		}	

		if ( y < m_height - 1 )
		{
			nPix = pix + m_width;
	
			if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
				e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
				                           m_smoothcost(beta_label,m_labeling[nPix]));

		}
	}
}


/**************************************************************************************/

void GCoptimization::setLabelOrder(bool RANDOM_LABEL_ORDER)
{
	m_random_label_order = RANDOM_LABEL_ORDER;
}
//<---Msun
void GCoptimization::setPairThingOrder(bool RANDOM_PAIRTHING_ORDER)
{
	m_random_pair_thing_order = RANDOM_PAIRTHING_ORDER;
}
// >---Msun

/****************************************************************************/
/* This procedure checks if an error has occured, terminates program if yes */

void GCoptimization::terminateOnError(bool error_condition,const char *message)

{ 
   if  (error_condition) 
   {
#ifdef MEX_COMPILE
       mexErrMsgIdAndTxt("GraphCut:internal_error","\nGCoptimization error: %s\n", message);
#else
       printf("GraphCut:internal_error","\nGCoptimization error: %s\n", message);
	   exit(1);
#endif
    }
}


