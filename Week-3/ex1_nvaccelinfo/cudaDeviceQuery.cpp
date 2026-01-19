#include <stdio.h>

int main(int argc, char *argv[])
{
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int device = 0; device < nDevices; device++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, device);
        int memClock;
        cudaDeviceGetAttribute(&memClock,cudaDevAttrMemoryClockRate, device);
        printf("Device %i: \"%s\".\n", device, prop.name);
        printf("  Multiprocessors: %i\n", prop.multiProcessorCount);
        printf("  Maximum number of threads per multiprocessor: %d.\n",
               prop.maxThreadsPerMultiProcessor);
        printf("  Peak Memory Bandwidth (GB/s): %f\n",
               2.0*memClock*(prop.memoryBusWidth/8)/1.0e6);
    }
}
