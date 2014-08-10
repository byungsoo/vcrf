/* GCoptimization.h */
/* Olga Veksler, 2005, olga@csd.uwo.ca */
/* Copyright 2005 Olga Veksler (olga@csd.uwo.ca)*/
/* email olga@csd.uwo.ca for any questions, suggestions and bug reports 

 
 /*****************    IMPORTANT!!!!!!!!!!!***********************************************************

  If you use this software, YOU HAVE TO REFERENCE (at least) 3 papers, the citations [1]
  [2] and [3] below

/****************************************************************************************************/

/******************************************************************************/
/*	For description and  example usage see GC_README.TXT                      */
/******************************************************************************/
/*

	This software library implements the Graph Cuts Energy Minimization methods
	described in 
	
	  
		[1] Efficient Approximate Energy Minimization via Graph Cuts 
		    Yuri Boykov, Olga Veksler, Ramin Zabih, 
		    IEEE transactions on PAMI, vol. 20, no. 12, p. 1222-1239, November 2001. 


    
	This software can be used only for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	If you wish to use this software (or the algorithms described in the aforementioned paper)
	for commercial purposes, you should be aware that there is a US patent: 

		R. Zabih, Y. Boykov, O. Veksler, 
		"System and method for fast approximate energy minimization via graph cuts ", 
		United Stated Patent 6,744,923, June 1, 2004



/* Together with the library implemented by O. Veksler, we provide, with the permission of the
    V. Kolmogorov and Y. Boykov the following libraries:  

  
1)	energy.h, which was developed by Vladimir Kolmogorov and  implements binary energy minimization 
	technique described in

	[2] What Energy Functions can be Minimized via Graph Cuts?
	    Vladimir Kolmogorov and Ramin Zabih. 
	    To appear in IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI). 
	    Earlier version appeared in European Conference on Computer Vision (ECCV), May 2002. 


	We use "energy.h" to implement the binary energy minimization step
	for the alpha-expansion and swap algorithms. The graph construction provided by "energy.h" is
	more efficient (and slightly more general) than the original graph construction for the 
	alpha-expansion algorithm in the paper cited as [1]
	

	This software can be used only for research purposes. IF YOU USE THIS SOFTWARE,
	YOU SHOULD CITE THE AFOREMENTIONED PAPER [2] IN ANY RESULTING PUBLICATION.

	

2) 	graph.h, block.h, maxflow.cpp
	
	This software library implements the maxflow algorithm
	described in

		[3] An Experimental Comparison of Min-Cut/Max-Flow Algorithms
		    for Energy Minimization in Vision.
		    Yuri Boykov and Vladimir Kolmogorov.
		    In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
		    September 2004

	This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
	at Siemens Corporate Research. To make it available for public use,
	it was later reimplemented by Vladimir Kolmogorov based on open publications.

	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication. 
*/

/* These 4 files (energy.h,graph.h, block.h, maxflow.cpp) are included in the curent library with permission of
Vladimir Kolmogorov and Yuri Boykov. They are included in a separate subdirectory named NotVekslerCode.
The can  also be downloaded independently from Vladimir Kolmogorov's
website:http://research.microsoft.com/~vnk/
Please Note that Vladimir's website is likely to move in the future                                         */


/****************************************************************************************************/
/*

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/



#ifndef __GCOPTIMIZATION_H__
#define __GCOPTIMIZATION_H__

#include "LinkedBlockList.h"
#include "bagon_assert.h" /* <assert.h>*/
#include "graph.h"
#include "energy.h"
#define m_datacost(pix,lab)     (m_datacost[(pix)*m_num_labels+(lab)] )
#define m_losscost(pix,lab)     (m_losscost[(pix)*m_num_labels+(lab)] ) //Msun 3/5 added loss cost
#define m_smoothcost(lab1,lab2) (m_smoothcost[(lab1)+(lab2)*m_num_labels] )
#define SET_INDIVIDUALLY 0
#define SET_ALL_AT_ONCE   1

// <!-- bagon 
#ifdef MEX_COMPILE
#include "mex.h" 
#include "GraphCut.h"
#endif
// bagon -->

// <!--Msun
#include "QPBO.h"
// Msun -->


class GCoptimization
{
public:


	///////////////////////////////////////////////////////////////////////////////////////
	//    First define needed data types                                                 //
	///////////////////////////////////////////////////////////////////////////////////////
	/* Type for the total energy calculation. Change it in Graph.h   */
	typedef Graph::flowtype EnergyType;

	/* Type for the individual terms in the energy function.Change it in Graph.h */
	typedef Graph::captype EnergyTermType;

	/* Type of label. Can be set to char, short, int, long */
	typedef int LabelType;

	/* Type for pixel. Can be set to  short, int, long */ 
	typedef int PixelType;


	/* The following 4 definitions are functions used to set data and smoothness terms */
	typedef EnergyTermType (*smoothFnPix)(PixelType pix1, PixelType pix2, LabelType lab1, LabelType lab2);
	typedef EnergyTermType (*smoothFnCoord)(PixelType x, PixelType y, LabelType l1, LabelType l2); 
	typedef EnergyTermType (*dataFnPix)(PixelType pix, LabelType lab);
	typedef EnergyTermType (*dataFnCoord)(PixelType x, PixelType y, LabelType l); 


	///////////////////////////////////////////////////////////////////////////////////////
	//    Next define constructors                                                       //
	///////////////////////////////////////////////////////////////////////////////////////


	/* Use this constructor only for grid of size width by height.  If you use this constructor,  */
	/* 4 connected grid neigbhorhood structure is assumed.  num_labels is the number of labels,   */
	/* Labels are assumed to be in  {0, 1, 2,....num_labels-1}                                    */
	/* dataSetup and smoothSetup can take only 2 values, namely SET_INDIVIDUALLY (0) and       */
	/* SET_ALL_AT_ONCE(1)                                                                       */ 
	/* In case dataSetup = SET_INDIVIDUALLY, to set the data cost you must use the             */
	/* member function setDataCost(pixel,label,cost).  If dataSetup = SET_ALL_AT_ONCE,          */
	/* to set the data costs you must use any of the setData() functions.                         */
	/* In case smoothSetup = SET_INDIVIDUALLY, to set the smooth cost you must use the         */
	/* member function setSmoothCost(label1,label2,cost).  If smoothSetup = SET_ALL_AT_ONCE,    */
	/* to set the smoothness costs you must use any of the setSmoothness() functions.             */
	GCoptimization(PixelType width,PixelType height,int num_labels,int dataSetup, int smoothSetup);


	/* This is the constructor for non-grid graphs. Everything is the same as in the constructor  */
	/* above, but neighborhood structure must  be specified by any of the two setNeighbors()      */
	/* functions */
	GCoptimization(PixelType num_pixels,int num_labels,int dataSetup, int smoothSetup);

	/* This constructor is the same as the first constructor, but array m_answer for                */
	/* storage of the labels is passed. This will save space, as the answer will not have to be     */
	/* stored twice, once in the current class, and once in the class which calls this optimization */
	/* class */
	GCoptimization(LabelType *m_answer,PixelType num_pixels,int nLabels,int dataSetup, int smoothSetup);

	/* This constructor is the same as the second constructor, but array m_answer for               */
	/* storage of the labels is passed. This will save space, as the answer will not have to be     */
	/* stored twice, once in the current class, and once in the class which calls this optimization */
	/* class */
	GCoptimization(LabelType *m_answer,PixelType width,PixelType height,int num_labels,int dataSetup, int smoothSetup);
	 
	// <!--Msun initialize the CRF with higher order terms
	void InitHOP(int nHigher);
	void InitDetHOP(int nHigher);
	void InitOccHOP(int nHigher);
	void InitPairThing(int nPairThing);
	void InitClassCo(bool useClassCoFlag);
	void setSegLoss(EnergyTermType* lossArray);
	void printThing();
	int setOneHOP(int n, int* ind, EnergyTermType* weights, EnergyTermType* gammas, EnergyTermType Q);
	//int setOneDemHOg(int n, int* ind, EnergyTermType* weights, EnergyTermType* gammas, EnergyTermType Q);
	int setOneDetHOP(int n, int* ind, EnergyTermType* weights, EnergyTermType* gammas, EnergyTermType Q, EnergyTermType thing_uw);
	int setOneOccHOP(int n, int* ind, EnergyTermType* weights, EnergyTermType* gammas, EnergyTermType Q, int ind2Det);
	void setOnePairThing( int pti, int* ind, EnergyTermType* weights);
	void setClassCo(  EnergyTermType * uw, EnergyTermType * pw);
	void setCompClassCo(  EnergyTermType * uw, EnergyTermType * pw);
	void InitCompLabel(int num_comp_labels);
	int setOneCompLabel(int n, int* ind);
	// Msun -->


	/* Peforms expansion algorithm. Runs the number of iterations specified by max_num_iterations */
	/* Returns the total energy of labeling   */
	EnergyType expansion(int max_num_iterations);

	/* Peforms expansion algorithm. Runs it until convergence */
	/* Returns the total energy of labeling   */
	EnergyType expansion();

	/* Peforms one iteration (one pass over all labels)  of expansion algorithm.*/
	EnergyType oneExpansionIteration();

	/* Peforms  expansion on one label, specified by the input parameter alpha_label */
	EnergyType alpha_expansion(LabelType alpha_label);

	/* Peforms  expansion on label alpha_label only for pixels specified by *pixels.  */
	/* num is the size of array pixels                                                */
	EnergyType alpha_expansion(LabelType alpha_label, PixelType *pixels, int num);
	
	/* Peforms swap algorithm. Runs it the specified number of iterations */
	EnergyType swap(int max_num_iterations);

	/* Peforms swap algorithm. Runs it until convergence */
	EnergyType swap();

	/* Peforms one iteration (one pass over all labels)  of swap algorithm.*/
	EnergyType oneSwapIteration();

	/* Peforms  swap on a pair of labels, specified by the input parameters alpha_label, beta_label */
	EnergyType alpha_beta_swap(LabelType alpha_label, LabelType beta_label);

	~GCoptimization();

	/* Returns current label assigned to input pixel */
	inline LabelType whatLabel(PixelType pixel){assert(pixel >= 0 && pixel < m_num_pixels);return(m_labeling[pixel]);};

	/* Sets data cost of pixel to label */
	inline void setDataCost(PixelType pixel, LabelType label, EnergyTermType cost){
		 assert(m_dataInput == SET_INDIVIDUALLY);
		 assert(pixel >= 0 && pixel < m_num_pixels && label >= 0 && label < m_num_labels );
	     m_datacost(pixel,label) = cost;};

	/* Sets smooth cost for label1, label2 to cost. */
	inline void setSmoothCost(LabelType label1, LabelType label2, EnergyTermType cost){
		assert(m_smoothInput == SET_INDIVIDUALLY);
		assert(label1 >= 0 && label1 < m_num_labels && label2 >= 0 && label2 < m_num_labels );
		m_smoothcost(label1,label2) = cost;};

	/* Makes pixel1 and pixel2 neighbors of each other. Can be called only 1 time for each         */
	/* unordered pair of pixels. Parameter weight can be used to set spacially varying terms       */
	/* If the desired penalty for neighboring pixels pixel1 and pixel2 is                          */
	/*  V(label1,label2) = weight*SmoothnessPenalty(label1,label2), then                           */
	/* member function setLabel should be called as: setLabel(pixel1,pixel2,weight)                */
	void setNeighbors(PixelType pixel1, PixelType pixel2, EnergyTermType weight);

	/* If the desired penalty for neighboring pixels pixel1 and pixel2 is                          */
	/*  V(label1,label2) = SmoothnessPenalty(label1,label2), then                                  */
	/* member function setLabel should be called as: setLabel(pixel1,pixel2)                       */
	/* Also use this function to set up the neighborhood structure when energy terms are specified */
	/* by a function pointer                                                                       */
	/* Again, this function can only be called one time for each distinct pair of pixels           */ 
	void setNeighbors(PixelType pixel1, PixelType pixel2);

	/* This function can be used to change the label of any pixel at any time      */
	inline void setLabel(PixelType pixel, LabelType label){
		assert(label >= 0 && label < m_num_labels && pixel >= 0 && pixel < m_num_pixels);m_labeling[pixel] = label;};
	
    // <!-- bagon
    /* Set all labels at once according to user's input                            */
    /* labels is an array of size num_pixels (or width*height in grid graph)       */
    void SetAllLabels(const LabelType* labels);
    
    /* copy the internal labels array into a given external array                  */
    void ExportLabels(LabelType* labels);
    
    /* get number of pixels in graph                                               */
    inline PixelType GetNumPixels() { return m_num_pixels; };
    
    /* get width, if not a grid returns one thus height = num_of_pixels / width */
    inline PixelType GetWidth() { 
        if (m_grid_graph) 
            return m_width; 
        return 1; };
    
    /* get number of possible labels                                               */
    inline LabelType GetNumLabels() { return m_num_labels; };
    // bagon -->
    
	/* Returns Data Energy of current labeling */
	EnergyType giveDataEnergy();

	/* Returns Smooth Energy of current labeling */
	EnergyType giveSmoothEnergy();

	//<--Msun
	/* Returns HOP Energy of current labeling */
	EnergyType compute_energy();
	EnergyType giveHOPEnergy();
	EnergyType giveClassDetHOP2Energy( LabelType class_label);
	EnergyType giveClassDetHOP2XYEnergy( LabelType class_label);
	EnergyType giveClassDetHOP2PriorEnergy( LabelType class_label);
	EnergyType giveDetHOP2Energy();
	EnergyType giveClassDetHOPEnergy( LabelType class_label);
	EnergyType giveDetHOPEnergy();
	EnergyType giveClassOccHOPEnergy( LabelType class_label);
	EnergyType giveOccHOPEnergy();
	EnergyType giveSpecPairThingEnergy( int thing1, int thing2);
	EnergyType givePairThingEnergy();
	EnergyType giveClassCoEnergy();
	EnergyType giveCompClassCoEnergy();
	void ExportDetLabels( LabelType* labels);
	void SetDetBoolLabels(const bool* det_bool_labels);
	void SetDetFixBool(const bool* det_bool_labels);
    EnergyType giveSegLossEnergyArray();
    inline int GetnHigherDet() { 
            return m_nHigherDet; 
    };

	void setPairThingOrder(bool RANDOM_DETBOOL_ORDER);
	//-->Msun

	/* By default, the labels are visited in random order for both the swap and alpha-expansion moves */
	/* Use this function with boolean argument 0 to fix the order to be not random                    */
	/* Use this function with argumnet 1 to fix the order back to random                              */
	void setLabelOrder(bool RANDOM_LABEL_ORDER);


	/* This function is used to set the data term, and it can be used only if dataSetup = SET_ALL_AT_ONCE */
	/* DataCost is an array s.t. the data cost for pixel p and  label l is stored at                        */
	/* DataCost[pixel*num_labels+l].  If the current neighborhood system is a grid, then                    */
	/* the data term for label l and pixel with coordinates (x,y) is stored at                              */ 
	/* DataCost[(x+y*width)*num_labels + l]. Thus the size of array DataCost is num_pixels*num_labels       */
	 void setData(EnergyTermType *DataCost);

	/* This function is used to set the data term, and it can be used only if dataSetup = SET_ALL_AT_ONCE */
	/* dataFn is a pointer to a function  f(Pixel p, Label l), s.t. the data cost for pixel p to have       */
	/* label l  is given by f(p,l) */
	 void setData(dataFnPix dataFn);


	 /* This function is used to set the data term, and it can be used only if dataSetup = SET_ALL_AT_ONCE */
	 /* and only for the grid graph                                                                          */
	 /* dataFn is a pointer to a function  f(x, y, l), s.t. the data cost for pixel with coordinates (x,y)   */
	 /* to have  label l  is given by f(x,y,l) */
	 void setData(dataFnCoord dataFn);


	 /* This function is used to set the smoothness term, and it can be used only if                         */
	 /* smoothSetup = SET_ALL_AT_ONCE                                                                      */
	 /*  V is an array of costs, such that V(label1,label2)  is stored at V[label1+num_labels*label2]        */
	 /* If graph is a grid, then using this  function only if the smooth costs are not spacially varying     */
	 /* that is the smoothness penalty V depends only on labels, but not on pixels                           */
	 void setSmoothness(EnergyTermType* V);

	 /* This function is used to set the smoothness term, and it can be used only if                         */
	 /* smoothSetup = SET_ALL_AT_ONCE AND the graph is a grid                                                */
     /* V is an array of costs, such that V(label1,label2)  is stored at V[label1+num_labels*label2]         */
	 /* hCue() is used to store the spacially varying coefficient w_pq.                                      */
	 /* if p has coordinates (x,y) and q has coordinates (x+1,y) then w_pq = hCue[x+y*width]                 */
	 /* The smoothness penalty, therefore, is                                                                */
	 /* V_pq(label1,label2) = V[label1+num_labels*label2]*hCue[x+y*width]                                    */
	 /* vCue() is used to store the spacially varying coefficient w_pq.                                      */
	 /* if p has coordinates (x,y) and q has coordinates (x,y+1) then w_pq = vCue[x+y*width]                 */
	 /* The smoothness penalty, therefore, is                                                                */
	 /* V_pq(label1,label2) = V[label1+num_labels*label2]*vCue[x+y*width]                                    */
	 void setSmoothness(EnergyTermType* V,EnergyTermType* hCue, EnergyTermType* vCue);

	 /* This function is used to set the smoothness term, and it can be used only if                         */
	 /* smoothSetup = SET_ALL_AT_ONCE                                                                      */
     /* cost is a function f(pix1,pix2,label1,label2) such that smoothness penalty for neigboring pixels     */
	 /* pix1 and pix2 to  have labels, respectively, label1 and label2 is f(pix1,pix2,label1,label2)         */
	 //void setSmoothness(smoothFnCoord cost);
	 void setSmoothness(smoothFnPix cost);

	 /* This function is used to set the smoothness term, and it can be used only if                         */
	 /* smoothSetup = SET_ALL_AT_ONCE AND the graph is a grid                                              */
     /* horz_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
	 /* (x,y) and (x+1,y) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
     /* vert_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
	 /* (x,y) and (x,y+1) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
	 void setSmoothness(smoothFnCoord horz_cost, smoothFnCoord vert_cost);

     // <!-- bagon
     /* validate class function */
#ifdef MEX_COMPILE
     bool IsClassValid() { return class_sig == VALID_CLASS_SIGNITURE; }; 
#else
     bool IsClassValid() { return class_sig == 0xabcd0123; }; 
#endif

     // bagon -->

	 // <!--Msun
	 void UseQPBO();
	 int m_num_labels;
	 PixelType m_num_pixels;
	 int m_num_comp_labels;
	 int m_hopType;
	 void SetSegWeight(EnergyTermType seg_weight);
	 void SetDetClassWeights(double* det_class_weights);
	 // Msun -->

private:

	//<<!--Msun quick hacks
	typedef enum 
	{
		GCsolver,
		QPBOsolver,
		NOIDEA,
	} solverType;
	solverType m_solver;
	// Msun -->

	typedef struct NeighborStruct{
		PixelType  to_node;
		EnergyTermType weight;
	} Neighbor;

	typedef enum 
	{
		ARRAY,
		FUNCTION_PIX,
		FUNCTION_COORD,
		NONE,
	} representationType;

    int class_sig; /* bagon: signiture value to verify class is ok */
    
	int m_width;
	int m_height;
	//PixelType m_num_pixels;
	LabelType *m_labeling;

	representationType m_dataType;
	representationType m_smoothType;
	bool m_random_label_order;
	bool m_grid_graph;
	bool m_varying_weights;               // this boolean flag is used only for grid graphs
	int  m_dataInput;
	int  m_smoothInput;
	bool m_deleteLabeling;

	EnergyTermType *m_datacost;
	EnergyTermType *m_losscost;//Msun: 3/5 added loss
	EnergyTermType *m_smoothcost;
	EnergyTermType *m_vertWeights;
	EnergyTermType *m_horizWeights;
	LinkedBlockList *m_neighbors;

	LabelType *m_labelTable;
	PixelType *m_lookupPixVar;
    
	EnergyTermType m_weight;

	/* Pointers to function for energy terms */
	dataFnPix m_dataFnPix;
	dataFnCoord m_dataFnCoord;

	smoothFnCoord m_horz_cost,m_vert_cost;
	smoothFnPix m_smoothFnPix;

	//<<!--Msun add higher-order-terms need for Eq.20 in pushmeet's IJCV paper
	PixelType *m_lookupValidPix;

	// --- regular_higher_order_potential---
	EnergyTermType *m_higherCost; // gamma_k values for the HOpotentials (#labels+1)x(#higher), last entry at each column is gamma_max
	EnergyTermType *m_higherTruncation; // Q_k values - one per HOpotential
	EnergyTermType **m_higherWeights; // BAGON: w_i weights for each node at a HOpotential
	EnergyTermType *m_higherP; // BAGON: sum w_i of each HOpotential (P value)
	int m_nHigher; // number of higher order potentials
	int **m_higherIndex; // an array per HOpotential with indices of participating nodes
 	int *m_higherElements; // for each HOpotential - the number of nodes participating in it, i.e., there are higherElements[j-1] nodes in HOpotential j

	// --- detection_higher_order_potential---
	bool *m_det_bool_label; // the current {0,1} label of the bool indicator
	bool *m_det_bool_negation; // true if construct the graph using negation of indicator
	PixelType *m_det_bool_lookupValid; // true if the indicator should be fixed inthe graph-cut, it is assgined -1, other is the variable indices used in the graph construction
	bool *m_det_fix_bool;// Msun added 9/13 to manually fix some detection indicators
	
	EnergyTermType *m_higherDetCost; // gamma_k values for the HOpotentials (#labels+1)x(#higher), last entry at each column is gamma_max
	LabelType *m_det_labeling;
	EnergyTermType *m_higherDetTruncation; // Q_k values - one per HOpotential
	EnergyTermType **m_higherDetWeights; // BAGON: w_i weights for each node at a HOpotential
	EnergyTermType *m_higherDetP; // BAGON: sum w_i of each HOpotential (P value)
	int m_nHigherDet; // number of higher order potentials
	int **m_higherDetIndex; // an array per HOpotential with indices of participating nodes
 	int *m_higherDetElements; // for each HOpotential - the number of nodes participating in it, i.e., there are higherElements[j-1] nodes in HOpotential j
	// Min 2/27 add new thing_unary_w
	EnergyTermType *m_higherDetUW;

	// --- occlusion_higher_order_potential---
	int *m_Occl2HigherDetIndex; // index for each occl potential corresponding to HigherDet
	EnergyTermType *m_higherOccCost; // gamma_k values for the HOpotentials (#labels+1)x(#higher), last entry at each column is gamma_max
	EnergyTermType *m_higherOccTruncation; // Q_k values - one per HOpotential
	EnergyTermType **m_higherOccWeights; // BAGON: w_i weights for each node at a HOpotential
	EnergyTermType *m_higherOccP; // BAGON: sum w_i of each HOpotential (P value)
	int m_nHigherOcc; // number of higher order potentials
	int **m_higherOccIndex; // an array per HOpotential with indices of participating nodes
 	int *m_higherOccElements; // for each HOpotential - the number of nodes participating in it, i.e., there are higherElements[j-1] nodes in HOpotential j

	// --- pair-thing potential ---
	int m_nPairThing; // number of Pair for things
	int **m_PairThingIndice; // [y_k y_l]
	EnergyTermType **m_PairThingWeights; //[w00 w01 w10 w11
	bool m_random_pair_thing_order;
	int *m_PairThingTable;

	// --- class co-occurrence potential ---
	bool m_useClassCoFlag; // if true should use class co-occurrence potential
	EnergyTermType *m_class_co_u_w; // class_co-occurrence unary weight
	EnergyTermType *m_class_co_p_w; // class_co-occurrence pair-wise weight
	bool *m_classExistFlag; // a cache of which classes exist in current m_labeling
	// Msun added 2/26
//	int m_num_comp_labels;
	LabelType *m_label2comp_label;
	int *m_num_dup_labels;
	LabelType **m_comp_label2label;
    // Msun added 3/5
    // 1) add learned wrights
    EnergyTermType m_seg_weight;
    EnergyTermType *m_det_class_weights;
	// Msun -->

	EnergyType start_expansion(int max_iterations);
	EnergyType start_swap(int max_iterations);

	void commonGridInitialization( PixelType width, PixelType height, int nLabels);
	void commonNonGridInitialization(PixelType num_pixels, int num_labels);
	void commonInitialization(int dataSetup, int smoothSetup);	

	void scramble_label_table();
	void perform_alpha_expansion(LabelType label);	
	// <!--Msun
	void scramble_pair_thing_table();
	bool CheckRegularity( int det_bool1, int det_bool2, int PtInd);
	PixelType SetDetBoolNegLookUp( LabelType alpha_label);
	void perform_alpha_expansion_QPBO(LabelType label);
	void set_up_expansion_energy_HOP( int size, LabelType alpha_label, Energy *e, Energy::Var *variables);
	void set_up_expansion_energy_Det_HOP2( LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables);
	void set_up_expansion_energy_Det_HOP( LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables);
	void set_up_expansion_energy_Occ_HOP( LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *det_variables);
	void set_up_expansion_energy_Pair_Thing( LabelType alpha_label, Energy *e, Energy::Var *det_variables, int det_size);
	void set_up_expansion_energy_Class_Co( LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *co_variables);
	void set_up_expansion_energy_CompClass_Co( LabelType alpha_label, Energy *e, Energy::Var *variables, Energy::Var *co_variables);
	void set_up_expansion_energy_QPBO_HOP( int size, LabelType alpha_label, QPBO<REAL> *e, int *variables);
	void set_up_expansion_energy_QPBO_DET( int size, LabelType alpha_label, QPBO<REAL> *e, int *variables);
	int getMaxLabel(int i);
	EnergyTermType cardinality(int i, int label);
	// Msun -->
	void perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label);

	EnergyType giveDataEnergyArray();
	EnergyType giveDataEnergyFnPix();
	EnergyType giveDataEnergyFnCoord();

	EnergyType giveSmoothEnergy_G_ARRAY_VW();
	EnergyType giveSmoothEnergy_G_ARRAY();
	EnergyType giveSmoothEnergy_G_FnPix();
	EnergyType giveSmoothEnergy_G_FnCoord();
	EnergyType giveSmoothEnergy_NG_ARRAY();
	EnergyType giveSmoothEnergy_NG_FnPix();

	void add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
	void add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
	void add_t_links_FnCoord(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
			
	// <!--Msun
	void add_t_links_ARRAY_QPBO(QPBO<REAL> *e, int *variables,int size,LabelType alpha_label);
	void add_t_links_FnPix_QPBO(QPBO<REAL> *e, int *variables,int size,LabelType alpha_label);
	void add_t_links_FnCoord_QPBO(QPBO<REAL> *e, int *variables,int size,LabelType alpha_label);
	// Msun -->

	void add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
	void add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
	void add_t_links_FnCoord_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
	
	void set_up_expansion_energy_G_ARRAY_VW(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
	void set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
	void set_up_expansion_energy_G_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
	void set_up_expansion_energy_G_FnCoord(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
	void set_up_expansion_energy_NG_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);		
	void set_up_expansion_energy_NG_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);		
	void set_up_expansion_energy_G_ARRAY_VW_pix(int size, LabelType alpha_label,Energy *e,
													  Energy::Var *variables, PixelType *pixels, int num );
	void set_up_expansion_energy_G_ARRAY_pix(int size, LabelType alpha_label,Energy *e,
													  Energy::Var *variables, PixelType *pixels, int num );

	// <!--Msun
	void set_up_expansion_energy_G_ARRAY_VW_QPBO(int size, LabelType alpha_label, QPBO<REAL>* e, int *variables);//
	void set_up_expansion_energy_G_ARRAY_QPBO(int size, LabelType alpha_label,QPBO<REAL>* e, int *variables);//
	void set_up_expansion_energy_G_FnPix_QPBO(int size, LabelType alpha_label,QPBO<REAL>* e, int *variables); //
	void set_up_expansion_energy_G_FnCoord_QPBO(int size, LabelType alpha_label,QPBO<REAL>* e, int *variables);//
	void set_up_expansion_energy_NG_ARRAY_QPBO(int size, LabelType alpha_label,QPBO<REAL>* e, int *variables);	//	
	void set_up_expansion_energy_NG_FnPix_QPBO(int size, LabelType alpha_label,QPBO<REAL>* e, int *variables);	//
	//void set_up_expansion_energy_G_ARRAY_VW_pix_QPBO(int size, LabelType alpha_label,Energy *e,
	//												  Energy::Var *variables, PixelType *pixels, int num );
	//void set_up_expansion_energy_G_ARRAY_pix_QPBO(int size, LabelType alpha_label,Energy *e,
	//												  Energy::Var *variables, PixelType *pixels, int num );
	int count_neighbor_Grid(int size, LabelType alpha_label);
	int count_neighbor_NGrid(int size, LabelType alpha_label);
	// Msun -->

	void set_up_swap_energy_G_ARRAY_VW(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
	void set_up_swap_energy_G_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
	void set_up_swap_energy_G_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
	void set_up_swap_energy_G_FnCoord(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
	void set_up_swap_energy_NG_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);		
	void set_up_swap_energy_NG_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);		


	void initialize_memory();
	void terminateOnError(bool error_condition,const char *message);

};

#endif

