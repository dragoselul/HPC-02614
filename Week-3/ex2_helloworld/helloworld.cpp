#include <omp.h>
#include <stdio.h>
int main(int argc, char *argv[]) 
{
    #pragma omp target teams parallel \
    num_teams(114) thread_limit(4*32)
    {
        int team_id = omp_get_team_num();
        int tid = omp_get_thread_num();

        int global_tid = team_id * 32 * 4 + tid;
        if (global_tid == 100) { int *a = (int*) 0x10000; *a = 0; }
        printf("Hello world from (%d, %d) global_tid: %d!\n",
        omp_get_team_num(),
        omp_get_thread_num(),global_tid);
    } // end target teams parallel
    return(0);
}