#include <iostream>
#include <omp.h>

using namespace std;


__global__ void kernel_call(int N, double *in, double* out)
{
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  
  for (int i = id; i < N; i += blockDim.x * gridDim.x)
    out[i] = in[i];
}


int main(){

    double *host_in, *host_out;
    double *dev_in, *dev_out;

    size_t N = 18874368; 
 
		
    //create buffer on host	
    host_in = (double*) malloc(N* sizeof(double));
    host_out = (double*) malloc(N * sizeof(double));


    //create buffer on device
    cudaError_t err = cudaMalloc(&dev_in, N*sizeof(double));
    if (err != cudaSuccess){
      cout<<"Dev Memory not allocated"<<endl;
      exit(-1);
    }


    err = cudaMalloc(&dev_out, N*sizeof(double));
    if (err != cudaSuccess){
       cout<<"Dev Memory not allocated"<<endl;
       exit(-1);
    }



    for (int i = 1; i < 128; i <<= 1)
{
 
   cout<<i<<" "<<N/i  <<" ";
    
    //using OpenMP to perform timing on the host   
    double st = omp_get_wtime();
    cudaMemcpy(dev_in, host_in, N * sizeof(double), cudaMemcpyHostToDevice);
    double et = omp_get_wtime();
 
    cout<<"Copy time: "<<(et-st)*1000<<"ms ";     


    //create GPU timing events for timing the GPU
    cudaEvent_t st2, et2;
    cudaEventCreate(&st2);
    cudaEventCreate(&et2);        
     
    cudaEventRecord(st2);
    kernel_call<<<1, 32>>>(N, dev_in, dev_out);
    cudaEventRecord(et2);
        
    //host waits until et2 has occured     
    cudaEventSynchronize(et2);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, st2, et2);

    cout<<"Kernel time: "<<milliseconds<<"ms"<<endl;

    //copy data out
    cudaMemcpy(host_out, dev_out, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaEventDestroy(st2);
    cudaEventDestroy(et2);
}
    free(host_in);
    free(host_out);
    cudaFree(dev_in);
    cudaFree(dev_out);

  return 0;
}
