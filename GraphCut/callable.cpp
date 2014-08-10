#include <iostream>
#include "callable.h"
#include "GCoptimization.h"
//#include "GraphCut.h"

#ifdef MEX_COMPILE
// something you will use in MEX
#else
// something you don't want to include in MEX complie eg. #include "mex.h"
#endif

using namespace std;


extern "C" {

void* GC_setLabelOrder(void* pInst) {
	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setLabelOrder( false);
}

void* GC_setPairThingOrder(void* pInst) {
	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setPairThingOrder( false);
}

void* GC_setSegLossCostMin(void* pInst, float* SegLoss) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setSegLoss(SegLoss);
}

void* GC_setDataCostMin(void* pInst, float* data) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	cout << "GC_setDataCostMin() is called." << endl;
	pInstance->setData( data);
	cout << "GC_setDataCostMin() is ended." << endl;
}

void* GC_setDataCostMin2(void* pInst, float* data, int num_pixels, int num_labels) {
	int i,j,c=0;

	float *data_array = new float[num_pixels * num_labels];

	for(j=0; j<num_pixels; j++) {
		for(i=0; i<num_labels; i++) {
			data_array[c] = data[c];
			c++;
		}
	}
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setData(data_array);

	delete [] data_array;
}
void* GC_setDataCost(void* pInst, double** data, int num_pixels, int num_labels) {
	int i,j,c;

//	cout << "GC_setData() is called." << endl;
	float *data_array = new float[num_pixels * num_labels];

	for(j=0; j<num_pixels; j++) {
		for(i=0; i<num_labels; i++) {
			//if(data[i][j] != 0) {
			//	cout << "data["<<i << "," <<j<<"]:"<<data[i][j]<<endl;
			//}
			
			data_array[j*num_labels + i] = (float)data[i][j];
//			cout << "data[i][j]:" << data[i][j] <<endl;
		}
	}
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setData(data_array);

	delete [] data_array;
}

void* GC_setSmoothCostMin(void* pInst, float *SmoothnessCost) {
	cout << "GC_setSmoothCostMin() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setSmoothness(SmoothnessCost);
}

void* GC_setSmoothCostMin2(void* pInst, float* smooth, int num_labels) {
	int i,j,c;

	cout << "GC_setSmoothness() is called." << endl;
	float *smooth_array = new float[num_labels * num_labels];
	c = 0;
	for(i=0; i<num_labels; i++) {
		for(j=0; j<num_labels; j++) {
			smooth_array[c] = smooth[c];
			c++;
		}
	}

  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setSmoothness(smooth_array);
	cout << "GC_setSmoothness() ended." << endl;

	delete [] smooth_array;
}

void* GC_setSmoothCost(void* pInst, double** smooth, int num_labels) {
	int i,j,c;

	cout << "GC_setSmoothness() is called." << endl;
	float *smooth_array = new float[num_labels * num_labels];
	c = 0;
	for(i=0; i<num_labels; i++) {
		for(j=0; j<num_labels; j++) {
			//cout << "smooth["<<i << "," <<j<<"]:"<<smooth[i][j]<<endl;
			smooth_array[c] = (float)smooth[i][j];
			c++;
		}
	}

  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setSmoothness(smooth_array);

	delete [] smooth_array;
}

void* GC_setNeighbors(void* pInst, int i, int j, double c) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->setNeighbors(i, j, (float)c);
}

void* GC_initHOP(void* pInst, int nHigher) {
	printf("GC_initHOP(%d) is called.\n", nHigher); fflush(stdout);
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->InitHOP(nHigher);
}

void* GC_initDetHOP(void* pInst, int nHigher) {
	printf("GC_initDetHOP(%d) is called.\n", nHigher); fflush(stdout);
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->InitDetHOP(nHigher);
}

void* GC_initOccHOP(void* pInst, int nHigher) {
	printf("GC_initOccHOP(%d) is called.\n", nHigher); fflush(stdout);
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->InitOccHOP(nHigher);
}

void* GC_initPairThing(void* pInst, int nHigher) {
	printf("GC_initPairThing(%d) is called.\n", nHigher); fflush(stdout);
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->InitPairThing(nHigher);
}

void* GC_initClassCo(void* pInst) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->InitClassCo(true);
}





void* GC_setLabels(void* pInst, int* labels) {
	printf("GC_setLabels() is called.\n"); fflush(stdout);
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->SetAllLabels(labels);
	printf("GC_setLabels() is done.\n"); fflush(stdout);

}


void* GC_setDetLabels(void* pInst, int* det_int_labels, int n_DetHOP) {
	
	printf("GC_setDetLabels(%d) is called.\n",n_DetHOP); fflush(stdout);

	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	bool* det_bool_labels;
	int i;

	det_bool_labels = new bool[n_DetHOP];
	

	for (i=0;i<n_DetHOP; i++) {
		if(det_int_labels[i] == 1) {
			det_bool_labels[i] = true;
		} else if(det_int_labels[i] == 0) {
			det_bool_labels[i] = false;
		} else {
			printf("det_int_labels[] should have integer 0 or 1 only.\n"); fflush(stdout);
			exit(1);
		}
	}
	printf("Int labels are converted into bool labels.\n"); fflush(stdout);
	pInstance->SetDetBoolLabels(det_bool_labels);
	delete [] det_bool_labels;

	printf("GC_setDetLabels() is done.\n"); fflush(stdout);
}



void* GC_setHOP(void* pInst, int nElement, int* ind, float* w, int num_gamma, float* gamma, float Q) {
	int *ind_new = new int[nElement];
	float *w_new = new float[nElement];
	float *gamma_new = new float[num_gamma];
	for (int ii=0; ii<nElement; ii++){
		ind_new[ii] = ind[ii]-1;// in c/c++ index start from 0
		w_new[ii] = w[ii];
	}
	for (int ii=0; ii<num_gamma; ii++){
		gamma_new[ii] = gamma[ii];
	}

  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	if(pInstance->setOneHOP(nElement, ind_new, w_new, gamma_new, Q)< 0)
		{printf("Failed to load hop.\n"); fflush(stdout);}

	delete ind_new;
	delete w_new;
	delete gamma_new;
}


void* GC_setDetHOP(void* pInst, int nElement, int* ind, float* w, int num_gamma, float* gamma, float Q, float uw) {
	int *ind_new = new int[nElement];
	float *w_new = new float[nElement];
	float *gamma_new = new float[num_gamma];
	for (int ii=0; ii<nElement; ii++){
		ind_new[ii] = ind[ii]-1; // in c/c++ index start from 0
		w_new[ii] = w[ii];
	}
	for (int ii=0; ii<num_gamma; ii++){
		gamma_new[ii] = gamma[ii];
	}

  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	if(pInstance->setOneDetHOP(nElement, ind_new, w_new, gamma_new, Q, uw)< 0)
		{printf("Failed to load hop.\n"); fflush(stdout);}
	
	delete ind_new;
	delete w_new;
	delete gamma_new;
}


void* GC_setOccHOP(void* pInst, int nElement, int* ind, float* w, float* gamma, float Q, int ind2Det) {
//	cout << "GC_setOccHOP() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	int i;
	if(pInstance->setOneOccHOP(nElement, ind, w, gamma, Q, ind2Det)< 0)
		{
			printf("Failed to load OccHOP., Q:%f\n", Q); fflush(stdout);
			
		}
}

void* GC_setPairThing(void* pInst, int idx, int* det_inds, float* w) {
//	cout << "GC_setPairThing() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	pInstance->setOnePairThing(idx, det_inds, w);
}

void* GC_initCompLabel(void *pInst, int numlabel) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	pInstance->InitCompLabel(numlabel);

}

void* GC_setOneCompModel(void *pInst, int n, int* ind) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	pInstance->setOneCompLabel(n, ind);

}


void* GC_setClassCo(void *pInst, double **approx_co_occur, int n_label, double weight) {

	int i, j, k;
	float* UnartClassCoCost = new float[n_label];
	float* PairClassCoCost = new float[n_label * n_label];
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	k = 0;
	for(i=0;i<n_label;i++) {
		UnartClassCoCost[i] = (float)weight * (float)approx_co_occur[i][i];
		for(j=0;j<n_label;j++) {
			PairClassCoCost[k] = (float)weight * (float)approx_co_occur[i][j];
			if(i==j) PairClassCoCost[k] = 0;
			k++;
		}
	}


	pInstance->setClassCo( UnartClassCoCost, PairClassCoCost);

	delete [] UnartClassCoCost;
	delete [] PairClassCoCost;
}

/* Msun: 3/5 added for getting Seg-based energy*/
double GC_compute_seg_energy(void* pInst) {

	GCoptimization::EnergyType e=0;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = pInstance->giveSmoothEnergy();
	e = e + pInstance->giveDataEnergy();
	e = e + pInstance->giveHOPEnergy();

	return e;
}

void GC_set_seg_weight( void* pInst, double weight){
	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->SetSegWeight( (GCoptimization::EnergyType) weight);
}

void GC_set_det_class_weights( void* pInst, double* weights){
	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->SetDetClassWeights( weights);
}
/* Msun: 3/5 added*/

double GC_compute_energy(void* pInst) {
//	cout << "GC_compute_energy() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = pInstance->compute_energy();;
	/*e = e + pInstance->giveSmoothEnergy();
	e = e + pInstance->giveDataEnergy();
	e = e + pInstance->giveHOPEnergy();

	printf("%f, %f, %f\n", pInstance->giveSmoothEnergy(),  pInstance->giveDataEnergy(),  pInstance->giveHOPEnergy());
	*/
	return e;
}

double GC_compute_energy_HOP(void* pInst) {
//	cout << "GC_compute_energy_HOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->giveHOPEnergy();

	return e;
}

double GC_compute_energy_DetHOP(void* pInst) {
//	cout << "GC_compute_energy_DetHOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	//e = e + pInstance->giveDetHOPEnergy();
	e = e + pInstance->giveDetHOP2Energy();

	return e;
}



double GC_compute_energy_ClassDetHOP(void* pInst, int class_label) {
//	cout << "GC_compute_energy_ClassDetHOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	//e = e + pInstance->giveClassDetHOPEnergy(class_label);
	e = e + pInstance->giveClassDetHOP2Energy(class_label);

	return e;
}


double GC_compute_energy_ClassOccHOP(void* pInst, int class_label) {
//	cout << "GC_compute_energy_ClassOccHOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->giveClassOccHOPEnergy(class_label);

	return e;
}


double GC_compute_energy_OccHOP(void* pInst) {
//	cout << "GC_compute_energy_OccHOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->giveOccHOPEnergy();
//	printf("GC_compute_energy_OccHOP() is done.\n"); fflush(stdout);

	return e;
}

double GC_compute_energy_PairThing(void* pInst) {
//	cout << "GC_compute_energy_PairThing() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->givePairThingEnergy();

	return e;
}

double GC_compute_energy_SpecPairThing(void *pInst, int thing1, int thing2) {
//	cout << "GC_compute_energy_PairThing() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->giveSpecPairThingEnergy(thing1, thing2);

	return e;
}

double GC_compute_energy_ClassCo(void* pInst) {
//	cout << "GC_compute_energy_HOP() is called." << endl;

	GCoptimization::EnergyType e;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	e = 0;
	e = e + pInstance->giveClassCoEnergy();

	return e;
}




void* GC_expansion(void* pInst) { //, int expansion_type) {
	//cout << "GC_expansion() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	pInstance->expansion(); //expansion_type);
}

void* GC_UseQPBO(void* pInst) {
	cout << "GC_UseQPBO() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	pInstance->UseQPBO();
	pInstance->m_hopType = 1;
}

void* GC_UseQPBO_DET(void* pInst) {
	cout << "GC_UseQPBO_DET() is called." << endl;
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);

	pInstance->UseQPBO();
	pInstance->m_hopType = 0;
}

int GC_whatLabel(void* pInst, int i) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	return pInstance->whatLabel(i);
}

void* GC_whatDetLabel(void* pInst, int* det_int_label, int n_DetHOP) {
  	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	//det_int_label = new int[n_DetHOP];

	pInstance->ExportDetLabels(det_int_label);
}

void* GCoptimizationWrapper(void *pInst)
{
  cout << "GCoptimizationWrapper() called" << endl;
  GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
  //pInstance->DoSomeAction();

  return NULL;
}

int GetNumPixels(void *pInst){
	cout<<"in GetNumPixels"<<endl;
	GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
	cout<<"in GetNumPixels done static_cast"<<endl;
	cout<<"GetNumPixels:m_num_pixels"<<pInstance->m_num_pixels<<endl;
	cout<<"GetNumPixels:m_num_labelss"<<pInstance->m_num_labels<<endl;
    return pInstance->m_num_pixels;
}

void* CreateGCoptimization(int num_sites,int num_labels) {
	cout << "CreateGCoptimization(n_site:"<<num_sites<<",n_labels:"<<num_labels<<") is called." << endl;
	GCoptimization* MyGraph = new GCoptimization(num_sites, num_labels, 1, 1);
	//cout<<"CreateGCoptimization:m_num_pixels"<<MyGraph->m_num_pixels<<endl;
	return static_cast<void*>(MyGraph);
}

void* DeleteGCoptimization(void *pInst) {
  //cout<<"in DeleteGCoptimization"<<endl;
  GCoptimization* pInstance = static_cast<GCoptimization*>(pInst);
  delete pInstance;
  pInstance = NULL;
  //cout<<"DeleteGCoptimization out"<<endl;
}


}
