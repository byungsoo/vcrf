
#include "mex.h"
#include "GCoptimization.h"
#include "GraphCut.h"
#include <stdlib.h>

/* Defines */


/*
 * Matlab wrapper for Weksler graph cut implementation
 *
 * usage:
 * [gch] = GraphCutConstrSparseHOP(dc, sc, SparseSc, hop)
 *
 * Note that data types are crucials!
 * 
 * Inputs:
 *  dc - of type float, array size [#labels*#nodes], the data term for node
 *             n recieving label l is stroed at [n*#labels + l]
 *  sc - of type float, array size [#labels.^2] the cost between l1 and l2 is
 *                   stored at [l1+l2*#labels] = Vpq(lp,lq)
 *  SparseSc - Sparse matrix of type double defining both the graph structure 
 *              and spatially dependent smoothness term
 *  hop - higher order potential array of structs with (#higher) entries, each entry:
 *     .ind - indices of nodes belonging to this hop
 *     .w - weights w_i for each participating node
 *     .gamma - #labels + 1 entries for gamma_1..gamma_max
 *     .Q - truncation value for this potential (assumes one Q for all labels)
 *  hop_det - higher order potential array of structs with (#higher) entries, each entry:
 *     .ind - indices of nodes belonging to this hop
 *     .w - weights w_i for each participating node
 *     .gamma - #labels + 1 entries for gamma_1..gamma_max
 *     .Q - truncation value for this potential (assumes one Q for all labels)
 *  hop_occ - higher order potential array of structs with (#higher) entries, each entry:
 *     .ind - indices of nodes belonging to this hop
 *     .w - weights w_i for each participating node
 *     .gamma - #labels + 1 entries for gamma_1..gamma_max
 *     .Q - truncation value for this potential (assumes one Q for all labels)
 *     .det_ind - index of hop_det indicator, refer to hop_det(det_ind)
 *  pair_thing - pair-wise potential for things
 *	   .det_inds - pair of indices of det_bool
 *     .w - [00 01 10 11]
 *  ClassCo_uw - unary weight for each class
 *  ClassCo_pw - pair-wise weights for all classes
 *
 * Outputs:
 *  gch - of type int32, graph cut handle - do NOT mess with it!
 */

/*const int HOP_N_OF_FIELDS(4); // expecting 4 fields for the HOpotentials struct
const char* HOP_FIELDS[HOP_N_OF_FIELDS] = {"ind", "w", "gamma", "Q"};
const int Occ_HOP_N_OF_FIELDS(5); // expecting 5 fields for the HOpotentials struct
const char* Occ_HOP_FIELDS[Occ_HOP_N_OF_FIELDS] = {"ind", "w", "gamma", "Q", "det_ind"};*/
template<class T>
void GetArr(const mxArray* x, T* arr, T bias = 0);
void SetHOPWrapper( GCoptimization* MyGraph, const mxArray *prhs[], int loc, int type);
void SetCompLabelWrapper( GCoptimization* MyGraph, const mxArray *prhs[], int loc);

void mexFunction(
    int		  nlhs, 	/* number of expected outputs */
    mxArray	  *plhs[],	/* mxArray output pointer array */
    int		  nrhs, 	/* number of inputs */
    const mxArray	  *prhs[]	/* mxArray input pointer array */
    )
{
#ifdef LOG_ON
	mexPrintf("nrhs=%d",nrhs);
#endif
    /* check number of arguments */
    if (nrhs != 10) {
        mexErrMsgIdAndTxt("GraphCut:NarginError","Wrong number of input argumnets");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("GraphCut:NargoutError","Wrong number of output argumnets");
    }
    /* check inputs */
    if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS ) {
        mexErrMsgIdAndTxt("GraphCut:DataCost", "DataCost argument is not of type float");
    }
    if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS ) {
        mexErrMsgIdAndTxt("GraphCut:SmoothnessCost", "SmoothnessCost argument is not of type float");
    }
    if (! mxIsSparse(prhs[2]) ) {
        mexErrMsgIdAndTxt("GraphCut:SmoothnessCost", "Graph Structure Matrix is not sparse");
    }
    
    GCoptimization::PixelType num_pixels;
    int num_labels;
    
    num_pixels = mxGetN(prhs[2]);
    if (mxGetM(prhs[2]) != num_pixels) {
        mexErrMsgIdAndTxt("GraphCut:SmoothnessCost", "Graph Structure Matrix is no square");
    }
    num_labels = mxGetNumberOfElements(prhs[0])/num_pixels;
    if (mxGetNumberOfElements(prhs[1]) != num_labels*num_labels) {
        mexErrMsgIdAndTxt("GraphCut:SmoothnessCost", "Size does not match number of labels");
    }
    
    /* construct the graph */
    GCoptimization* MyGraph = new GCoptimization(num_pixels, num_labels, SET_ALL_AT_ONCE, SET_ALL_AT_ONCE);
    
    /* set the nieghbors and weights according to sparse matrix */
    
    double   *pr;
    mwIndex  *ir, *jc;
    mwSize   col, total=0;
    mwIndex  starting_row_index, stopping_row_index, current_row_index;

  
    /* Get the starting positions of all four data arrays. */
    pr = mxGetPr(prhs[2]); // Msun: prhs[2] is SparseSc (symetric sparse matrix)
    ir = mxGetIr(prhs[2]);
    jc = mxGetJc(prhs[2]);
    
    for (col=0; col<num_pixels; col++)  {
        starting_row_index = jc[col];
        stopping_row_index = jc[col+1];
        if (starting_row_index == stopping_row_index) {
            continue;
        } else {
            for (current_row_index = starting_row_index;
                current_row_index < stopping_row_index;
                current_row_index++)  {
                    /* use only upper triangle of matrix */
                    if ( ir[current_row_index] >= col ) {
                        MyGraph->setNeighbors(ir[current_row_index], col, pr[total++]);
                    } else {
                        total++;
                    }
                    
            }
        }
    }

	Graph::captype *DataCost = (Graph::captype*)mxGetData(prhs[0]);
    Graph::captype *SmoothnessCost = (Graph::captype*)mxGetData(prhs[1]);  
	 
    /* set data term */
    MyGraph->setData(DataCost);
    /* set the smoothness term */
    MyGraph->setSmoothness(SmoothnessCost);
    
	/* set the general hop data ================================*/
	SetHOPWrapper( MyGraph, prhs, 3, 0);
	/* END of set the general hop data ================================*/
    
	/* set the detection hop data ================================*/
	if (nrhs >4){
#ifdef LOG_ON
		mexPrintf("Det HOP\n");
#endif
		SetHOPWrapper( MyGraph, prhs, 4, 2);
	}
	/* END of set the detection hop data ================================*/

	/* set the detection hop data ================================*/
	if (nrhs >5){
#ifdef LOG_ON
		mexPrintf("Occ HOP is not functioning now\n");
#endif
//		SetHOPWrapper( MyGraph, prhs, 5, 2);
	}
	/* END of set the detection hop data ================================*/

	/* set pair-wise thing potential =============================*/
	if (nrhs >6){
		const int N_OF_FIELDS = 2;
		const char* FIELDS[2] = {"det_inds", "w"};
		int nPairThing,ii;
		int fields_indices[N_OF_FIELDS]; // indices to fields in hop struct
		if ( mxGetClassID(prhs[6]) != mxSTRUCT_CLASS )
			mexErrMsgIdAndTxt("robustpn:inputs","hop must be a struct array");
		nPairThing = mxGetNumberOfElements(prhs[6]);
#ifdef LOG_ON
		mexPrintf("nPairThing=%d\n",nPairThing);
#endif
		if (nPairThing!=0 ){
			MyGraph->InitPairThing(nPairThing);
			// expecting FIELDS fieds
			if ( mxGetNumberOfFields(prhs[6]) != N_OF_FIELDS )
				mexErrMsgIdAndTxt("robustpn:inputs","hop must have %d fields", N_OF_FIELDS);
			// chack that we have the right fields
			for ( ii = 0; ii < N_OF_FIELDS ; ii++ ) {
				fields_indices[ii] = mxGetFieldNumber(prhs[6], FIELDS[ii]);
				if ( fields_indices[ii] < 0 )
					mexErrMsgIdAndTxt("robustpn:inputs","hop is missing %s field", FIELDS[ii]);
			}
			// Add the HO potentials
			mxArray *xdet_inds, *xw;
			int * det_inds, n;
			Graph::captype* w;
			for ( ii = 0 ; ii < nPairThing; ii++ ) {
				xdet_inds = mxGetFieldByNumber(prhs[6], ii, fields_indices[0]);
				n = mxGetNumberOfElements(xdet_inds);
				//mexPrintf("num_ele=%d\n",n);
				det_inds = new int[n]; // allocation for energy
				GetArr(xdet_inds, det_inds, -1); // bias = -1 convert from 1-ind of matlab to 0-ind of C
				
				xw = mxGetFieldByNumber(prhs[6], ii, fields_indices[1]);
				n = mxGetNumberOfElements(xw);
				w = new Graph::captype[n]; // allocation for energy
				GetArr(xw, w);

				MyGraph->setOnePairThing(ii, det_inds, w);
			}
		}
	}
	/* END of set pair-wise thing potential =============================*/

	/* set clas co-occurrence potential =============================*/	
	if (nrhs==9){
		int num = mxGetNumberOfElements(prhs[7]);
		//if (num == MyGraph->m_num_labels){

		if (num != 0){
#ifdef LOG_ON
			mexPrintf("use Class Co");
#endif
			Graph::captype *UnartClassCoCost = (Graph::captype*)mxGetData(prhs[7]);
			Graph::captype *PairClassCoCost = (Graph::captype*)mxGetData(prhs[8]);  
			MyGraph->InitClassCo(true);
			MyGraph->setClassCo( UnartClassCoCost, PairClassCoCost);
		}
	}	
	/* END set clas co-occurrence potential =============================*/	
	if (nrhs>9) {
		int num = mxGetNumberOfElements(prhs[7]);

		if (num != 0){
#ifdef LOG_ON
			mexPrintf("use CompClass Co");
#endif
			Graph::captype *UnartClassCoCost = (Graph::captype*)mxGetData(prhs[7]);
			Graph::captype *PairClassCoCost = (Graph::captype*)mxGetData(prhs[8]);  
			MyGraph->InitClassCo(true);
			SetCompLabelWrapper( MyGraph, prhs, 9);
			MyGraph->setCompClassCo( UnartClassCoCost, PairClassCoCost);
		}
	}

	/* set Compress Label used by co-occurrence potential =============================*/	

	/* END set Compress Label used by co-occurrence potential =============================*/	

#ifdef LOG_ON
    mexPrintf("Done reading hops\n");
#endif
        
    /* create a container for the pointer */
    const mwSize dims[2] = {1,0};
    plhs[0] = mxCreateNumericArray(1, /*(int*)*/dims, MATLAB_POINTER_TYPE, mxREAL);
    
    GraphHandle* gh;
    gh = (GraphHandle*) mxGetData(plhs[0]);
    *gh = (GraphHandle)(MyGraph);
}

template<class T>
void GetArr(const mxArray* x, T* arr, T bias)
{
    int ii, n = mxGetNumberOfElements(x);
    void *p = mxGetData(x);
    char* cp;
    unsigned char* ucp;
    short* sp;
    unsigned short* usp;
    int* ip;
    unsigned int* uip;
    int64_T *i64p;
    uint64_T *ui64p;
    double* dp;
    float* fp;
    
    switch (mxGetClassID(x)) {
        case mxCHAR_CLASS:
        case mxINT8_CLASS:    
            cp = (char*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = cp[ii] + bias;
            return;
        case mxDOUBLE_CLASS:
            dp = (double*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = dp[ii]+ bias;
            return;
        case mxSINGLE_CLASS:
            fp = (float*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = fp[ii]+ bias;
            return;
        case mxUINT8_CLASS:
            ucp = (unsigned char*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ucp[ii]+ bias;
            return;
        case mxINT16_CLASS:
            sp = (short*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = sp[ii]+ bias;
            return;
        case mxUINT16_CLASS:
            usp = (unsigned short*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = usp[ii]+ bias;
            return;
        case mxINT32_CLASS:
            ip = (int*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ip[ii]+ bias;
            return;
        case mxUINT32_CLASS:
            uip = (unsigned int*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = uip[ii]+ bias;
            return;
        case mxINT64_CLASS:
            i64p = (int64_T*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = i64p[ii]+ bias;
            return;
        case mxUINT64_CLASS:
            ui64p = (uint64_T*)p;
            for ( ii = 0 ; ii < n ; ii++)
                arr[ii] = ui64p[ii]+ bias;
            return;
        default:
            mexErrMsgIdAndTxt("robustpn:GetArr","unsupported data type");
    }
}
void SetHOPWrapper( GCoptimization* MyGraph, const mxArray *prhs[], int loc, int type){

	int N_OF_FIELDS;
	if (type == 0 || type == 1){
		N_OF_FIELDS = 4;
	}else{
		N_OF_FIELDS = 5;
	}
	//const char* FIELDS[5] = {"ind", "w", "gamma", "Q", "det_ind"}; for OccDet. obsolete
	const char* FIELDS[5] = {"ind", "w", "gamma", "Q", "thing_uw"};
	
	int nHigher, ii;
    int hop_fields_indices[N_OF_FIELDS]; // indices to fields in hop struct
    // check pin[3] is struct array with proper feilds
    if ( mxGetClassID(prhs[loc]) != mxSTRUCT_CLASS )
        mexErrMsgIdAndTxt("robustpn:inputs","hop must be a struct array");
	nHigher = mxGetNumberOfElements(prhs[loc]);
#ifdef LOG_ON
	mexPrintf("nHigher=%d\n",nHigher);
#endif
	if (nHigher ==0 ){
		return;
	}
	if (type == 0){
		MyGraph->InitHOP(nHigher);
	}else if (type == 1){
		MyGraph->InitDetHOP(nHigher);
	}else{
		//MyGraph->InitOccHOP(nHigher); for OccDet. obsolete
		MyGraph->InitDetHOP(nHigher);
	}
    // expecting FIELDS fieds
    if ( mxGetNumberOfFields(prhs[loc]) != N_OF_FIELDS )
        mexErrMsgIdAndTxt("robustpn:inputs","hop must have %d fields", N_OF_FIELDS);
    // chack that we have the right fields
    for ( ii = 0; ii < N_OF_FIELDS ; ii++ ) {
        hop_fields_indices[ii] = mxGetFieldNumber(prhs[loc], FIELDS[ii]);
        if ( hop_fields_indices[ii] < 0 )
            mexErrMsgIdAndTxt("robustpn:inputs","hop is missing %s field", FIELDS[ii]);
    }
    // Add the HO potentials
    //mxArray *xind, *xw, *xgamma, *xQ, *xdet_ind; for OccDet. obsolete
    mxArray *xind, *xw, *xgamma, *xQ, *x_thing_uw;
    int * ind, n;
    Graph::captype* w;
    Graph::captype* gamma;
    Graph::captype Q;
	//int det_ind; for OccDet. obsolete
	Graph::captype thing_uw;
    for ( ii = 0 ; ii < nHigher; ii++ ) {
        xind = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[0]);
        n = mxGetNumberOfElements(xind);
		//mexPrintf("num_ele=%d\n",n);
        ind = new int[n]; // allocation for energy
        GetArr(xind, ind, -1); // bias = -1 convert from 1-ind of matlab to 0-ind of C
        
        xw = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[1]);
        if ( mxGetNumberOfElements(xw) != n ) {
#ifdef LOG_ON
			mexPrintf("n=%d w_length=%d\n",n,mxGetNumberOfElements(xw));
#endif
            delete[] ind;
            mexErrMsgIdAndTxt("robustpn:inputs","hop %d: number of indices is different than number of weights", ii);
        }
        w = new Graph::captype[n]; // allocation for energy
        GetArr(xw, w);
        
        xgamma = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[2]);
        if ( mxGetNumberOfElements(xgamma) != MyGraph->m_num_labels+1 ) {
            delete[] ind;
            delete[] w;
            mexErrMsgIdAndTxt("robustpn:inputs","hop %d: must have exactly %d gamma values", ii, MyGraph->m_num_labels+1);
        }
        gamma = new Graph::captype[ MyGraph->m_num_labels+1];
        GetArr(xgamma, gamma);
        
        xQ = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[3]);
        Q = (Graph::captype)mxGetScalar(xQ);

		if (N_OF_FIELDS > 4){
			x_thing_uw = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[4]);
			if ( mxGetNumberOfElements(x_thing_uw) != 1 ) {
#ifdef LOG_ON
				mexPrintf("n=%d w_length=%d\n",1,mxGetNumberOfElements(x_thing_uw));
#endif
            	delete[] ind;
            	delete[] w;
				delete[] gamma;
				mexErrMsgIdAndTxt("robustpn:inputs","hop %d: number of indices is different than number of weights", ii);
			}
			thing_uw = (Graph::captype)mxGetScalar(x_thing_uw);
		// for OccDet. obsolete
        //	xdet_ind = mxGetFieldByNumber(prhs[loc], ii, hop_fields_indices[4]);
        //	det_ind = (int)mxGetScalar(xdet_ind)-1; // bais= -1 covert from 1-ind of matlab to 0-ind of C
		}

		if (type == 0){
			if (MyGraph->setOneHOP(n, ind, w, gamma, Q) < 0 ) {
				delete[] gamma;
				mexErrMsgIdAndTxt("robustpn:inputs","failed to load hop #%d", ii);
			}
		}else if (type == 1){
//			if (MyGraph->setOneDetHOP(n, ind, w, gamma, Q) < 0 ) {
//				delete[] gamma;
//				mexErrMsgIdAndTxt("robustpn:inputs","failed to load hop #%d", ii);
//			}
		}else{
			//if (MyGraph->setOneOccHOP(n, ind, w, gamma, Q, det_ind) < 0 ) {
			if (MyGraph->setOneDetHOP(n, ind, w, gamma, Q, thing_uw) < 0 ) {
				delete[] gamma;
				mexErrMsgIdAndTxt("robustpn:inputs","failed to load hop #%d", ii);
			}
		}
        delete[] gamma; // this array is being allocated inside energy
        //mexPrintf("Done reading hop(%d) / %d\n", ii, nHigher);
    }

}

void SetCompLabelWrapper( GCoptimization* MyGraph, const mxArray *prhs[], int loc){

	int N_OF_FIELDS=1;
	const char* FIELDS[1] = {"ind"};
	
	int num_comp_labels, ii;
    int comp_label_fields_indices[N_OF_FIELDS]; // indices to fields in comp_label struct
    // check pin[3] is struct array with proper feilds
    if ( mxGetClassID(prhs[loc]) != mxSTRUCT_CLASS )
        mexErrMsgIdAndTxt("robustpn:inputs","comp_label must be a struct array");
	num_comp_labels = mxGetNumberOfElements(prhs[loc]);
#ifdef LOG_ON
	mexPrintf("num_comp_labels=%d\n",num_comp_labels);
#endif
	if (num_comp_labels ==0 ){
		return;
	}
	MyGraph->InitCompLabel(num_comp_labels);
    // expecting FIELDS fieds
    if ( mxGetNumberOfFields(prhs[loc]) != N_OF_FIELDS )
        mexErrMsgIdAndTxt("robustpn:inputs","comp_label must have %d fields", N_OF_FIELDS);
    // chack that we have the right fields
    for ( ii = 0; ii < N_OF_FIELDS ; ii++ ) {
        comp_label_fields_indices[ii] = mxGetFieldNumber(prhs[loc], FIELDS[ii]);
        if ( comp_label_fields_indices[ii] < 0 )
            mexErrMsgIdAndTxt("robustpn:inputs","comp_label is missing %s field", FIELDS[ii]);
    }
    // Add the HO potentials
    mxArray *xind;
    int * ind, n;
    for ( ii = 0 ; ii < num_comp_labels; ii++ ) {
        xind = mxGetFieldByNumber(prhs[loc], ii, comp_label_fields_indices[0]);
        n = mxGetNumberOfElements(xind);
		//mexPrintf("num_ele=%d\n",n);
        ind = new int[n]; // allocation for energy
        GetArr(xind, ind, -1); // bias = -1 convert from 1-ind of matlab to 0-ind of C
        
		if (MyGraph->setOneCompLabel(n, ind) < 0 ) {
			mexErrMsgIdAndTxt("robustpn:inputs","failed to load comp_label #%d", ii);
		}
    }
}
