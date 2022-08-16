/*
 * MEX-CUDA GPU Parallel Computing Code for Long Short Term Memory Forward Process
 * This uses the MATLAB CUDA API in a MEX function that takes gpuArray inputs and 
 * returns gpuArray outputs.
 *
 * Written by Sungwook Wi, UMass
 */


#include "mex.h"
#include "gpu/mxGPUArray.h"

/*
 * Device code
 */
 
void __global__ lstm_forward_hidden(float * const Hall, float * const Call,
									float * const Gi, float * const Gf, float * const Go, float * const Gg,
									float const * const Ui, float const * const Uf, float const * const Uo, float const * const Ug,
									float const * const Wi, float const * const Wf, float const * const Wo, float const * const Wg,
									float const * const bi, float const * const bf, float const * const bo, float const * const bg,
									float const * const X, float * const H, float * const C,
									int const NumIn, int const NumHid, int const NumSeq, int const BatchSize)
{
	
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int co_x, co_y, j, k, N, M;

	N = NumHid * BatchSize;	
	M = NumIn * BatchSize;
	
	for (k = 0; k < NumSeq; k++) {		
		
		if (i < N) {	
			
			co_x = i % NumHid;
			co_y = i / NumHid;
			Gi[i] = 0.0;
			Gf[i] = 0.0;
			Go[i] = 0.0;
			Gg[i] = 0.0;
			for (j=0; j < NumIn; j++) {
				Gi[i] += Ui[NumHid*j+co_x] * X[NumIn*co_y+j+k*M];
				Gf[i] += Uf[NumHid*j+co_x] * X[NumIn*co_y+j+k*M];
				Go[i] += Uo[NumHid*j+co_x] * X[NumIn*co_y+j+k*M];
				Gg[i] += Ug[NumHid*j+co_x] * X[NumIn*co_y+j+k*M];
			}

			co_x = i % NumHid;
			co_y = i / NumHid;
			for (j=0; j < NumHid; j++) {
				Gi[i] += Wi[NumHid*j+co_x] * H[NumHid*co_y+j];
				Gf[i] += Wf[NumHid*j+co_x] * H[NumHid*co_y+j];
				Go[i] += Wo[NumHid*j+co_x] * H[NumHid*co_y+j];
				Gg[i] += Wg[NumHid*j+co_x] * H[NumHid*co_y+j];						
			}
			
			Gi[i] += bi[co_x];
			Gf[i] += bf[co_x];
			Go[i] += bo[co_x];
			Gg[i] += bg[co_x];
			
			Gi[i] = 1.0 / (1.0 + expf(-Gi[i]));
			Gf[i] = 1.0 / (1.0 + expf(-Gf[i]));
			Go[i] = 1.0 / (1.0 + expf(-Go[i]));
			Gg[i] = tanhf(Gg[i]);
			

			C[i] = Gg[i]*Gi[i] + Gf[i]*C[i]; 			
			H[i] = tanhf(C[i])*Go[i];
			
			Call[i+N*k] = C[i];
			Hall[i+N*k] = H[i];
			
		}  
		
	}		

	
}	

void __global__ lstm_forward_out(float * const Yhat, float const * const V, float const * const bv, float * const H, 
									int const NumOut, int const NumHid, int const BatchSize,
									float const * const drop_rate, float const * const drop_ind)
{
	
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int co_x, co_y, j, N;
	float drop_scale_factor;
	
	drop_scale_factor = 1-drop_rate[0];

	N = NumOut * BatchSize;
	
	if (i < N) {	 

		co_x = i % NumOut;
		co_y = i / NumOut;
		Yhat[i] = 0.0;
		for (j=0; j < NumHid; j++) {
			if (drop_ind[j] != 1) { 
				Yhat[i] += V[NumOut*j+co_x] * H[NumHid*co_y+j]/drop_scale_factor;	
			}	
			__syncthreads();	
		}
		Yhat[i] += bv[co_x];
		
	}																				
							
}

void __global__ lstm_forward_loss(float * const Yhat, float const * const Yobs, float const * const err_scale, float * const Loss, 
									int const NumOut, int const BatchSize)
{
	
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int N;

	N = NumOut * BatchSize;
	
	if (i < N) {	 
		Loss[i] = 0.5*pow(Yhat[i]-Yobs[i],2)*err_scale[i];
	}																				
	
}
 

/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare all variables.*/
    // INPUTS
	mxGPUArray const *Ui;
	mxGPUArray const *Uf;
	mxGPUArray const *Uo;
	mxGPUArray const *Ug;
	mxGPUArray const *Wi;
	mxGPUArray const *Wf;
	mxGPUArray const *Wo;
	mxGPUArray const *Wg;
	mxGPUArray const *bi;
	mxGPUArray const *bf;
	mxGPUArray const *bo;
	mxGPUArray const *bg;
	mxGPUArray const *V;
	mxGPUArray const *bv;
	mxGPUArray const *X;
	mxGPUArray const *H;
	mxGPUArray const *C;		
	mxGPUArray const *Yobs;
	mxGPUArray const *err_scale;
	mxGPUArray const *dropRate;
	mxGPUArray const *dropInd;
    float const *d_Ui;
	float const *d_Uf;
	float const *d_Uo;
	float const *d_Ug;
	float const *d_Wi;
	float const *d_Wf;
	float const *d_Wo;
	float const *d_Wg;
	float const *d_bi;
	float const *d_bf;
	float const *d_bo;
	float const *d_bg;
	float const *d_V;
	float const *d_bv;
	float const *d_X;
	float *d_H;
	float *d_C;	
	float const *d_Yobs;
	float const *d_err_scale;
	float const *d_dropRate;
	float const *d_dropInd;
	
	// OUTPUTS
	mxGPUArray *Gi;
	mxGPUArray *Gf;
	mxGPUArray *Go;
	mxGPUArray *Gg;
	mxGPUArray *Call;
	mxGPUArray *Hall;
	mxGPUArray *Yhat;
	mxGPUArray *Loss;
	float *d_Gi;
	float *d_Gf;
	float *d_Go;
	float *d_Gg;
	float *d_Call;
	float *d_Hall;
	float *d_Yhat;
	float *d_Loss;

	
    int N;
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg1 = "Invalid input to MEX file: 20 inputs must be provided.";
    char const * const errMsg2 = "Invalid input to MEX file: Input(s) is not a GPU array.";
    char const * const errMsg3 = "Invalid input to MEX file: Input(s) is not a float array.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 1024;
    int blocksPerGrid;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
	if ((nrhs!=21)) {
        mexErrMsgIdAndTxt(errId, errMsg1);
    }

    if (!(mxIsGPUArray(prhs[0])) || 
		!(mxIsGPUArray(prhs[1])) ||
		!(mxIsGPUArray(prhs[2])) ||
		!(mxIsGPUArray(prhs[3])) ||
		!(mxIsGPUArray(prhs[4])) ||
		!(mxIsGPUArray(prhs[5])) ||
		!(mxIsGPUArray(prhs[6])) ||
		!(mxIsGPUArray(prhs[7])) ||
		!(mxIsGPUArray(prhs[8])) ||
		!(mxIsGPUArray(prhs[9])) ||
		!(mxIsGPUArray(prhs[10])) ||
		!(mxIsGPUArray(prhs[11])) ||
		!(mxIsGPUArray(prhs[12])) ||
		!(mxIsGPUArray(prhs[13])) ||
		!(mxIsGPUArray(prhs[14])) ||
		!(mxIsGPUArray(prhs[15])) ||
		!(mxIsGPUArray(prhs[16])) ||
		!(mxIsGPUArray(prhs[17])) ||
		!(mxIsGPUArray(prhs[18])) ||
		!(mxIsGPUArray(prhs[19])) ||
		!(mxIsGPUArray(prhs[20]))) 
    {
        mexErrMsgIdAndTxt(errId, errMsg2);
    }

	
    Ui = mxGPUCreateFromMxArray(prhs[0]);
	Uf = mxGPUCreateFromMxArray(prhs[1]);
	Uo = mxGPUCreateFromMxArray(prhs[2]);
	Ug = mxGPUCreateFromMxArray(prhs[3]);
	Wi = mxGPUCreateFromMxArray(prhs[4]);
	Wf = mxGPUCreateFromMxArray(prhs[5]);
	Wo = mxGPUCreateFromMxArray(prhs[6]);
	Wg = mxGPUCreateFromMxArray(prhs[7]);
	bi = mxGPUCreateFromMxArray(prhs[8]);
	bf = mxGPUCreateFromMxArray(prhs[9]);
	bo = mxGPUCreateFromMxArray(prhs[10]);
	bg = mxGPUCreateFromMxArray(prhs[11]);
	V = mxGPUCreateFromMxArray(prhs[12]);
	bv = mxGPUCreateFromMxArray(prhs[13]);
	X = mxGPUCreateFromMxArray(prhs[14]);	
	H = mxGPUCreateFromMxArray(prhs[15]);
	C = mxGPUCreateFromMxArray(prhs[16]);	
	Yobs = mxGPUCreateFromMxArray(prhs[17]);
	err_scale = mxGPUCreateFromMxArray(prhs[18]);
	dropRate = mxGPUCreateFromMxArray(prhs[19]);
	dropInd = mxGPUCreateFromMxArray(prhs[20]);
	

    /*
     * Verify that input really is a float array before extracting the pointer.
     */
    if ((mxGPUGetClassID(Ui) != mxSINGLE_CLASS) || 
    	(mxGPUGetClassID(Uf) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Uo) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Ug) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Wi) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Wf) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Wo) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Wg) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(bi) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(bf) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(bo) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(bg) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(V)  != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(bv) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(X)  != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(H)  != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(C)  != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(Yobs) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(err_scale) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(dropRate) != mxSINGLE_CLASS) ||
    	(mxGPUGetClassID(dropInd) != mxSINGLE_CLASS)) 
    {
        mexErrMsgIdAndTxt(errId, errMsg3);
    }
	
	

    /*
     * Now that we have verified the data type, extract a pointer to the input data on the device.
     */
	d_Ui = (float const *)(mxGPUGetDataReadOnly(Ui));
	d_Uf = (float const *)(mxGPUGetDataReadOnly(Uf));
	d_Uo = (float const *)(mxGPUGetDataReadOnly(Uo));
	d_Ug = (float const *)(mxGPUGetDataReadOnly(Ug));
	d_Wi = (float const *)(mxGPUGetDataReadOnly(Wi));
	d_Wf = (float const *)(mxGPUGetDataReadOnly(Wf));
	d_Wo = (float const *)(mxGPUGetDataReadOnly(Wo));
	d_Wg = (float const *)(mxGPUGetDataReadOnly(Wg));
	d_bi = (float const *)(mxGPUGetDataReadOnly(bi));
	d_bf = (float const *)(mxGPUGetDataReadOnly(bf));
	d_bo = (float const *)(mxGPUGetDataReadOnly(bo));
	d_bg = (float const *)(mxGPUGetDataReadOnly(bg));
	d_V = (float const *)(mxGPUGetDataReadOnly(V));
	d_bv = (float const *)(mxGPUGetDataReadOnly(bv));
	d_X = (float const *)(mxGPUGetDataReadOnly(X));	
	d_H = (float *)(mxGPUGetDataReadOnly(H));
	d_C = (float *)(mxGPUGetDataReadOnly(C));	
	d_Yobs = (float const *)(mxGPUGetDataReadOnly(Yobs));
	d_err_scale = (float const *)(mxGPUGetDataReadOnly(err_scale));
	d_dropRate = (float const *)(mxGPUGetDataReadOnly(dropRate));
	d_dropInd = (float const *)(mxGPUGetDataReadOnly(dropInd));


    /* Create a GPUArray to hold the result and get its underlying pointer. */
	const mwSize *dimsU = mxGPUGetDimensions(Ui);
	const mwSize *dimsX = mxGPUGetDimensions(X);	
	const mwSize *dimsV = mxGPUGetDimensions(V);
	size_t NumHid = dimsU[0];    // Number of hidden units
	size_t NumIn = dimsU[1];     // Number of input units (i.e., number of features)	
	size_t NumOut = dimsV[0];    // Number of output units	
	size_t NumSeq = dimsX[2];    // Sequence length
	size_t BatchSize = dimsX[1]; // Batch size (i.e., number of samples in a mini-batch)
	
	mwSize dims1[2] = {NumOut, BatchSize};
	mwSize dims2[2] = {NumHid, BatchSize};	
	mwSize dims3[3] = {NumHid, BatchSize, NumSeq};	

    Gi = mxGPUCreateGPUArray(2,dims2,mxGPUGetClassID(Ui),mxGPUGetComplexity(Ui),MX_GPU_DO_NOT_INITIALIZE);
	Gf = mxGPUCreateGPUArray(2,dims2,mxGPUGetClassID(Uf),mxGPUGetComplexity(Uf),MX_GPU_DO_NOT_INITIALIZE);
	Go = mxGPUCreateGPUArray(2,dims2,mxGPUGetClassID(Uo),mxGPUGetComplexity(Uo),MX_GPU_DO_NOT_INITIALIZE);
	Gg = mxGPUCreateGPUArray(2,dims2,mxGPUGetClassID(Ug),mxGPUGetComplexity(Ug),MX_GPU_DO_NOT_INITIALIZE);
	Call = mxGPUCreateGPUArray(3,dims3,mxGPUGetClassID(Ug),mxGPUGetComplexity(Ug),MX_GPU_DO_NOT_INITIALIZE);
	Hall = mxGPUCreateGPUArray(3,dims3,mxGPUGetClassID(Ug),mxGPUGetComplexity(Ug),MX_GPU_DO_NOT_INITIALIZE);
	Yhat = mxGPUCreateGPUArray(2,dims1,mxGPUGetClassID(Ug),mxGPUGetComplexity(Ug),MX_GPU_DO_NOT_INITIALIZE);
	Loss = mxGPUCreateGPUArray(2,dims1,mxGPUGetClassID(Ug),mxGPUGetComplexity(Ug),MX_GPU_DO_NOT_INITIALIZE);
    d_Gi = (float *)(mxGPUGetData(Gi));
	d_Gf = (float *)(mxGPUGetData(Gf));
	d_Go = (float *)(mxGPUGetData(Go));
	d_Gg = (float *)(mxGPUGetData(Gg));
	d_Call = (float *)(mxGPUGetData(Call));
	d_Hall = (float *)(mxGPUGetData(Hall));
	d_Yhat = (float *)(mxGPUGetData(Yhat));
	d_Loss = (float *)(mxGPUGetData(Loss));


    /*
     * Call the kernel using the CUDA runtime API. 
     * Here, I am using 1-d for grid and block deimension configuration.
     */

	N = (int)(NumHid * BatchSize);
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	lstm_forward_hidden<<<blocksPerGrid, threadsPerBlock>>>(d_Hall, d_Call, 
															d_Gi, d_Gf, d_Go, d_Gg, 
															d_Ui, d_Uf, d_Uo, d_Ug, 
															d_Wi, d_Wf, d_Wo, d_Wg,
															d_bi, d_bf, d_bo, d_bg,
															d_X, d_H, d_C, 
															NumIn, NumHid, NumSeq, BatchSize);
					
	N = (int)(NumOut * BatchSize);
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;					
	lstm_forward_out<<<blocksPerGrid, threadsPerBlock>>>(d_Yhat, d_V, d_bv, d_H, NumOut, NumHid, BatchSize, d_dropRate, d_dropInd);	
	
	N = (int)(NumOut * BatchSize);
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	lstm_forward_loss<<<blocksPerGrid, threadsPerBlock>>>(d_Yhat, d_Yobs, d_err_scale, d_Loss, NumOut, BatchSize);
	
	

    /* Wrap the result up as a MATLAB gpuArray for return. */  
	plhs[0] = mxGPUCreateMxArrayOnGPU(Yhat);
	plhs[1] = mxGPUCreateMxArrayOnGPU(Loss);
	plhs[2] = mxGPUCreateMxArrayOnGPU(Hall);
	plhs[3] = mxGPUCreateMxArrayOnGPU(Call);	
	
	

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    // Freeing INPUT Array
	mxGPUDestroyGPUArray(Ui);
	mxGPUDestroyGPUArray(Uf);
	mxGPUDestroyGPUArray(Uo);
	mxGPUDestroyGPUArray(Ug);
	mxGPUDestroyGPUArray(Wi);
	mxGPUDestroyGPUArray(Wf);
	mxGPUDestroyGPUArray(Wo);
	mxGPUDestroyGPUArray(Wg);
	mxGPUDestroyGPUArray(bi);
	mxGPUDestroyGPUArray(bf);
	mxGPUDestroyGPUArray(bo);
	mxGPUDestroyGPUArray(bg);
    mxGPUDestroyGPUArray(X);	
    mxGPUDestroyGPUArray(H);	
	mxGPUDestroyGPUArray(C);
	mxGPUDestroyGPUArray(V);	
	mxGPUDestroyGPUArray(bv);
	mxGPUDestroyGPUArray(Yobs);
	mxGPUDestroyGPUArray(err_scale);
	mxGPUDestroyGPUArray(dropRate);
	mxGPUDestroyGPUArray(dropInd);
	
	// Freeing OUTPUY Array
	mxGPUDestroyGPUArray(Gi);
	mxGPUDestroyGPUArray(Gf);
	mxGPUDestroyGPUArray(Go);
	mxGPUDestroyGPUArray(Gg);
	mxGPUDestroyGPUArray(Call);
	mxGPUDestroyGPUArray(Hall);
	mxGPUDestroyGPUArray(Yhat);
	mxGPUDestroyGPUArray(Loss);

}
