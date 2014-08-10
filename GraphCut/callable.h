#ifndef __CALLABLE_H
#define __CALLABLE_H

#ifdef __cplusplus
extern "C" {
#endif
	void* GC_setLabelOrder(void *pInst);
	void* GC_setPairThingOrder(void* pInst);
	void* GCoptimizationWrapper(void *pInst);
	void* GC_setSegLossCost(void* pInst, double** SegLoss, int num_pixels, int num_labels);
	void* GC_setDataCost(void* pInst, double** data, int num_pixels, int num_labels);
	void* GC_setDataCostMin2(void* pInst, float* data, int num_pixels, int num_labels);
	void* GC_setDataCostMin(void* pInst, float* data);
	void* GC_setSmoothCost(void* pInst, double** smooth, int num_labels);
	void* GC_setSmoothCostMin2(void* pInst, float* smooth, int num_labels);
	void* GC_setSmoothCostMin(void* pInst, float *SmoothnessCost);
	void* GC_setNeighbors(void* pInst, int i, int j, double c);
	void* GC_setLabels(void* pInst, int *labels);
	void* GC_setDetLabels(void* pInst, int* det_int_labels, int n_DetHOP);
	void* GC_initHOP(void* pInst, int nHigher);
	void* GC_initDetHOP(void* pInst, int nHigher);
	void* GC_initOccHOP(void* pInst, int nHigher);
	void* GC_initPairThing(void* pInst, int nHigher);
	void* GC_initClassCo(void* pInst);
	void* GC_initCompLabel(void *pInst, int numlabel);
	void* GC_setOneCompModel(void *pInst, int n, int* ind);

	void* GC_setHOP(void* pInst, int nElement, int* ind, float* w, int num_gamma, float* gamma, float Q);
	void* GC_setDetHOP(void* pInst, int nElement, int* ind, float* w, int num_gamma, float* gamma, float Q, float uw);
	void* GC_setOccHOP(void* pInst, int nElement, int* ind, float* w, float* gamma, float Q, int ind2Det);
	void* GC_setPairThing(void* pInst, int idx, int* det_inds, float* w);
	void* GC_setClassCo(void *pInst, double **approx_co_occur, int n_label, double weight);

	void* GC_UseQPBO(void* pInst);
	void* GC_UseQPBO_DET(void* pInst);

	double GC_compute_energy_HOP(void* pInst);
	double GC_compute_energy_DetHOP(void* pInst);
	double GC_compute_energy_ClassDetHOP(void* pInst, int class_label);
	double GC_compute_energy_ClassOccHOP(void* pInst, int class_label);
	double GC_compute_energy_OccHOP(void* pInst);
	double GC_compute_energy_PairThing(void* pInst);
	double GC_compute_energy_SpecPairThing(void* pInst, int thing1, int thing2);
	double GC_compute_energy_ClassCo(void* pInst);
	double GC_compute_energy(void* pInst);
    double GC_compute_seg_energy(void* pInst);
    void GC_set_seg_weight( void* pInst, double weight);
    void GC_set_det_class_weights( void* pInst, double* det_class_weights);
	//void* GC_expansion(void* pInst, int expansion_type);
	void* GC_expansion(void* pInst); //, int expansion_type);
	int GC_whatLabel(void* pInst, int i);
	void* GC_whatDetLabel(void* pInst, int* det_int_label, int n_DetHOP);
	int GetNumPixels(void *pInst);
	void* CreateGCoptimization(int num_sites,int num_labels);
	void* DeleteGCoptimization(void *pInst);
#ifdef __cplusplus
}
#endif

#endif

