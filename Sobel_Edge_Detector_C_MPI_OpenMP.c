#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

//----------------------------------------------------------time function-------------------------------------------
/*   ttype: type to use for representing time */
typedef double ttype;
ttype tdiff(struct timespec a, struct timespec b)
/* Find the time difference. */
{
    ttype dt = ((b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) / 1E9);
    return dt;
}

struct timespec now()
    /* Return the current time. */
{
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return t;
}

//clock_t begin, end;
struct timespec begin, end;
double time_spent;
//-----------------------------------------------------------------------------------------------------------------

//------------------------------------------declare parameters---------------------------------------------------------
int N;
int rank, size, tid, nthreads;
int rem;
int* rec_buf, * sendcounts, * displs, * update_rec_buf, * Ap;
int global_max, global_min;


//**************************************************************************************************************


//---------------------------------------------main function----------------------------------------------------
//void initialize_data(int A[N][N], int N);

void distribute_data(int A[N][N], int N);
void mask_operation(int A[N][N], int N);
void collect_results( int N);

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);                     //MPI init
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //------------------------------------------parameters---------------------------------------------------------
    // N = atof(argv[1]);
    int N = 5000;
    int A[N][N];    // origianl array


    //**************************************************************************************************************

    /*******begin read image data from txt to array*/
    FILE* readInput;
    readInput = fopen("5000testinput.txt", "r");
    //readInput = fopen("50testinput.txt", "r");
    if (readInput == NULL) {
        printf("File could not opened!");
    }

    int i, j, num;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            fscanf(readInput, "%d", &A[i][j]);

        }
    }
    fclose(readInput);
    /********************end reading data to array****************/

    begin = now();
    distribute_data (A, N); // use scatterv
    //test printf rank and scatterv output
    /*for (int i = 0; i < size; i++)
    {
        printf("%d.sendcount = %d \n", i, sendcounts[i]);
        printf("%d.displs = %d \n", i, displs[i]);
    }
    printf("rank : %d \n", rank);*/
    /*int count = 0;
    for( int i = 0; i < N; i++)
    {
        for( int j = 0; j < N; j++)
        {

            printf("scatter output: %d \t", rec_buf[i][j]);
            count++;
            if(count == sendcounts[rank])
            {break;}
        }

    }*/
    mask_operation(A, N);
    collect_results(N); // use gatherv
    //free(sendcounts);       // if you free memory others will get interrupt
    //free(displs);


    /****normalized threshold and denoise****/
    if (rank ==0)
    {
      printf("use MPI+omp rank %d: the global max is %d, the global min is %d\n", rank, global_max, global_min);
      float ratio = 0.0;
      /****add dynamic schedule for threshold and denoise****/
      /*int chunksize3 = 100;
      printf("chunksize3 is %d\n", chunksize3);*/

      #pragma omp parallel shared(nthreads, N, Ap) private(i,j,tid)
      {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        if(tid == 0)
        printf("rank: %d num of threads = %d\n", rank, nthreads);

        #pragma omp for
        //#pragma omp for schedule(dynamic, chunksize3)
        for(i = 1; i < N-1; i++) {
          for(j = 1; j < N-1; j++) {
            ratio = (double) (Ap[i*N+j] - global_min) / (global_max - global_min);
            Ap[i*N+j] = ratio * 255;
            if(Ap[i*N+j]<51)
            {
              Ap[i*N+j] = 0;
            }
          }
        }
      }
    }
    /****end normalized threshold and denoise****/
    end = now();
    time_spent = tdiff(begin, end);
    printf("rank %d whole process spend: %lf\n",rank, time_spent);
    //--------------------------write arr 2 txt------------------------
    if (rank == 0)
    {

        FILE* writeOutput = fopen("Part2out_MPIreduce.txt", "wb");

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if (j < N - 1)
                {
                    fprintf(writeOutput, "%d\t", Ap[i * N + j]);
                }
                else {
                    fprintf(writeOutput, "%d\t\n", Ap[i * N + j]);
                }
            }
        }
        fclose(writeOutput);
    }
    //------------------------------------------------------------------
    MPI_Finalize();
    return 0;

}

//**************************************************************************************************************

//------------------------------------------functions---------------------------------------------------------

void distribute_data(int A[N][N], int N)
{
        int i, j;

        sendcounts = (int*)malloc(sizeof(int) * size);
        // array describing how many elements to send to each process
        displs = (int*)malloc(sizeof(int) * size);
        // array describing the displacements where each segment begins
        rec_buf = (int*)malloc(sizeof(int) * (N * N));
        //init sendcounts, displs

        // calculate sendcounts and displs
        int row_proc, rem_row_proc, sum;
        row_proc = (N - 2) / size;
        rem_row_proc = (N - 2) % size;

        if ((N - 2) <= size)
        {
            for (i = 0; i < size; i++)
            {
                sendcounts[i] = 3 * N;
                displs[i] = i * N;
                if (i >= (N - 2))
                {
                    sendcounts[i] = 0;
                    displs[i] = 0;
                }
            }
        }

        if ((N - 2) > size)
        {
            sum = 0;
            for (i = 0; i < size; i++)
            {

                if (rem_row_proc == 0)
                {

                    sendcounts[i] = (row_proc + 2) * N;
                    displs[i] = sum * N;
                    sum = sum + row_proc;
                }
                if (rem_row_proc > 0)
                {
                    displs[i] = sum * N;
                    sendcounts[i] = (row_proc + 3) * N;
                    sum = sum + row_proc + 1;
                    rem_row_proc--;
                }
            }
        }printf("\n");


        printf("\n"); // print sendcount and displs
        if (rank == 0)
        {
            printf("scatter's sendcounts and displacement\n");
            for (i = 0; i < size; i++)
            {

                printf("sendcount[%d] = %d\n", i, sendcounts[i]);
                printf("displs[%d] = %d\n", i, displs[i]);
            }
        }


      MPI_Scatterv(&A[0][0], sendcounts, displs, MPI_INT, rec_buf, N*N, MPI_INT, 0, MPI_COMM_WORLD);
      printf("\n");

          /*printf("# %d rank receive matrix after scatterv(partial)\n", rank);
          printf("\n");
          for (i = 0; i < 15; i++)
          {
              for (j = 0; j < 15; j++)
              {
                  printf("%d\t", rec_buf[i * 15 + j]);
              }
              printf("\n");
          }*/



}

void mask_operation(int A[N][N], int N)
{
    int i, j;
    int qx, qy;
    update_rec_buf = (int*)malloc(sizeof(int) * (N * N));
    printf("\n");


    printf("rank id: %d \n", rank);
    // arrange rec_buf location and calculate value

    int chunksize1 = 100;
    if(rank ==0){
            printf("chunksize1 is %d\n", chunksize1);
            }
#pragma omp parallel shared(rank, update_rec_buf,chunksize1) private(tid,i,j, qx, qy)
    {

            tid = omp_get_thread_num();
            nthreads = omp_get_num_threads();
            if(tid == 0)
            printf("rank: %d num of threads = %d\n", rank, nthreads);

       int mx[3][3] = {
	    	 {-1, 0, 1},
		     {-2, 0, 2},
	   	   {-1, 0, 1}
        };
       int my[3][3] = {
	     	{1, 2, 1},
	    	{0, 0, 0},
	    	{-1, -2, -1}
      	};

        #pragma omp for schedule(dynamic, chunksize1)
        //#pragma omp for schedule(static,chunksize1)
        //#pragma omp for
        for (i = 0; i < sendcounts[rank] / N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if (rank == 0)
                {
                    if (i != 0 && i != ((sendcounts[rank] / N) - 1) && j != 0 && j != N - 1)
                    {

                        qx = ((-rec_buf[(i - 1) * N + j - 1] - rec_buf[(i + 1) * N + j - 1] - 2 * rec_buf[(i) * N + j -1]) +(2*rec_buf[(i) * N + j+1]
                            + rec_buf[(i - 1) * N + j + 1] + rec_buf[(i + 1) * N + j + 1]));
                        qy = ((rec_buf[(i - 1) * N + j - 1] + rec_buf[(i - 1) * N + j + 1] + 2 * rec_buf[(i-1) * N + j]) -2*rec_buf[(i+1) * N + j]
                            - rec_buf[(i + 1) * N + j - 1] - rec_buf[(i + 1) * N + j + 1]);
                        update_rec_buf[i * N + j] = abs(qx) + abs(qy);

                    }
                    else
                    {
                        update_rec_buf[i * N + j] = rec_buf[i * N + j];
                    }
                }
                if (rank == size - 1)
                {
                    if (i != 0 && i != ((sendcounts[rank] / N) - 1) && j != 0 && j != N - 1)
                    {

                        qx = ((-rec_buf[(i - 1) * N + j - 1] - rec_buf[(i + 1) * N + j - 1] - 2 * rec_buf[(i) * N + j -1]) +(2*rec_buf[(i) * N + j+1]
                            + rec_buf[(i - 1) * N + j + 1] + rec_buf[(i + 1) * N + j + 1]));
                        qy = ((rec_buf[(i - 1) * N + j - 1] + rec_buf[(i - 1) * N + j + 1] + 2 * rec_buf[(i-1) * N + j]) -2*rec_buf[(i+1) * N + j]
                            - rec_buf[(i + 1) * N + j - 1] - rec_buf[(i + 1) * N + j + 1]);
                        update_rec_buf[(i - 1) * N + j] = abs(qx) + abs(qy);
                    }

                    if ((i == (sendcounts[rank] / N - 2)) || (j == 0) || (j == N - 1))
                    {
                        update_rec_buf[i * N + j] = rec_buf[(i + 1) * N + j];
                    }


                }
                if (rank != 0 && rank != size - 1)
                {
                    if (i != 0 && i != ((sendcounts[rank] / N) - 1) && j != 0 && j != N - 1)
                    {

                        qx = ((-rec_buf[(i - 1) * N + j - 1] - rec_buf[(i + 1) * N + j - 1] - 2 * rec_buf[(i) * N + j -1]) +(2*rec_buf[(i) * N + j+1]
                            + rec_buf[(i - 1) * N + j + 1] + rec_buf[(i + 1) * N + j + 1]));
                        qy = ((rec_buf[(i - 1) * N + j - 1] + rec_buf[(i - 1) * N + j + 1] + 2 * rec_buf[(i-1) * N + j]) -2*rec_buf[(i+1) * N + j]
                            - rec_buf[(i + 1) * N + j - 1] - rec_buf[(i + 1) * N + j + 1]);
                        update_rec_buf[(i - 1) * N + j] = abs(qx) + abs(qy);
                    }
                    else
                    {
                        update_rec_buf[i * N + j] = rec_buf[(i + 1) * N + j];
                    }
                }

            }
        }
        printf("(calculate matrix)thread id = %d\n", tid);
    }
    /*                    verify it works
    if (rank == 0)
    {
        printf("\n");
        printf("original matrix in distribute function before scatterv \n");
        for (i = 0; i < N; i++)
        {//test
            for (j = 0; j < N; j++)
            {
                printf("%d \t", A[i][j]);
            }
            printf("\n");
        }
    }
            printf("\n");

            printf("this is from mask operation from rank %d\n",rank);
            for (i = 0; i < (sendcounts[rank] / 10); i++)
            {
                for (j = 0; j < 10; j++)
                {
                    printf("%d\t", update_rec_buf[i*10+j]);
                }
                printf("\n");
            }
     */
}

void collect_results(int N)
{
    int i, j;
    Ap = (int*)malloc(sizeof(int) * (N * N));


    //recalculate sendcounts and displs
    for (i = 0; i < size; i++)
    {
        if (i == 0)
        {

            sendcounts[i] = sendcounts[i] - N;

        }
        if (i == (size - 1))
        {
            sendcounts[i] = sendcounts[i] - N;
            displs[i] = displs[i] + N;
        }
        if (i != 0 && i != (size - 1))
        {
            sendcounts[i] = sendcounts[i] - 2 * N;
            displs[i] = displs[i] + N;
        }
    }

    /**********normalized threshold every rank works**********/
    int min = 5000, max = 0;  //reduction(max:max) reduction(min:min)
    /**********add dynamic schedule**********/
    /*int chunksize2 = 100;
    if(rank ==0){
            printf("chunksize2 is %d\n", chunksize2);
          }*/

    #pragma omp parallel shared(rank,N) private(tid, i, j)
    {
      if(rank == 0)  //ignore first line and edge column
      {
        #pragma omp for reduction(max:max) reduction(min:min)
        //#pragma omp for schedule(dynamic, chunksize2) reduction(max:max) reduction(min:min)
        for(i = 1; i < (sendcounts[rank]/N); i++)
        {
          for(j = 1; j < N-1; j++) {
            if (update_rec_buf[i*N+j] < min) {
              min = update_rec_buf[i*N+j];
            }
            else if (update_rec_buf[i*N+j] > max) {
              max = update_rec_buf[i*N+j];
            }
          }
        }
      }
      if(rank == size -1)   //ignore last line and edge column
      {
        #pragma omp for reduction(max:max) reduction(min:min)
        //#pragma omp for schedule(dynamic, chunksize2) reduction(max:max) reduction(min:min)
        for(i = 0; i < ((sendcounts[rank]/N)-1); i++)
        {
          for(j = 1; j < N-1; j++) {
            if (update_rec_buf[i*N+j] < min) {
              min = update_rec_buf[i*N+j];
            }
            else if (update_rec_buf[i*N+j] > max) {
              max = update_rec_buf[i*N+j];
            }
          }
        }
      }
      if (rank != 0 && rank != size - 1)   //ignore edge column
      {
        #pragma omp for reduction(max:max) reduction(min:min)
        //#pragma omp for schedule(dynamic, chunksize2) reduction(max:max) reduction(min:min)
        for(i = 0; i < (sendcounts[rank]/N); i++)
        {
          for(j = 1; j < N-1; j++) {
            if (update_rec_buf[i*N+j] < min) {
              min = update_rec_buf[i*N+j];
            }
            else if (update_rec_buf[i*N+j] > max) {
              max = update_rec_buf[i*N+j];
            }
          }
        }
      }
    }

    if (rank == 0)
    {
        printf("\n");
        printf("gather's sendcounts and displacements\n");
        for (i = 0; i < size; i++)
        {

            printf("sendcount[%d] = %d\n", i, sendcounts[i]);
            printf("displs[%d] = %d\n", i, displs[i]);
        }
    }

    MPI_Gatherv(update_rec_buf, sendcounts[rank], MPI_INT, Ap, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    /****use MPI_Reduce pass local max and min to master rank***/
    global_max = 0;
    global_min = 0;
    MPI_Reduce(&max,&global_max,1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&min,&global_min,1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    /****end use MPI_Reduce pass local max and min to master rank***/

    /*         test gather function
    if (rank == 0)
    {

        printf("\n");
        printf("gatherv receive(partial)\n");
        for (i = 0; i < 15; i++)
        {
            for (j = 0; j < 15; j++)
            {
                printf("%d\t", Ap[i * 15 + j]);
            }
            printf("\n");
        }
    }
    */
}
