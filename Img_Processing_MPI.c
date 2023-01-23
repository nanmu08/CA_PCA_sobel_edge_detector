/*collective_MPI.c nanmu integer type with all print and time print*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef double ttype;
ttype tdiff(struct timespec a, struct timespec b)
/* Find the time difference. */
{
  ttype dt = (( b.tv_sec - a.tv_sec ) + ( b.tv_nsec - a.tv_nsec ) / 1E9);
  return dt;
}

struct timespec now()
/* Return the current time. */
{
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t;
}

#define MASTER 0 //define rank 0 as master rank
int *Ap;
int taskid, size;
int i, j, send_element, calc_row, extra_row, max_element;
int *send_row, *sendcounts, *displs;
int *rec_buf; //for distribute func
int *result; //for mask_operation func
int *send_num, *send_buf; //for collect_results func
struct timespec begin, end;
double time_spent;
void initialize_data (int *A, int N);
void distribute_data (int *A, int N);
void mask_operation (int *A, int N, int*Ap);
void collect_results (int *Ap, int N);


/***********begin main function**********************************/
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int N = atof(argv[1]);
  int A[N*N]; // add matrix A (if use A[N*N], then function use func(A, N)), use A[N][N] func(A[0],N)
  initialize_data(A, N);
  distribute_data (A, N); // use scatterv
  mask_operation(A, N, Ap);
  collect_results(Ap, N); // use gatherv
  printf("-------------------------------------------------------\n\n");
  printf("rank %d total time is %.8f sec\n", taskid, time_spent);
  printf("-------------------------------------------------------\n\n");
  MPI_Finalize();
  return 0;
}

/**********************begin four functions*******************************/

void initialize_data (int *A, int N)
{
  //begin = clock();
  begin = now();
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(taskid == MASTER)
  {
    int seed = 1;
    srand(seed);
    //printf("master rank %d initialize matrix A\n", taskid);
    for (i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
      {
        A[i*N+j] = rand() % (256 - 0) +0;
        //printf("%d\t", A[i*N+j]);
      }

      //printf("\n");
    }
  }
  //MPI_Barrier (MPI_COMM_WORLD);
}

void distribute_data (int *A, int N)
{
  
  int begin_segm = 0;
  // array describing how many elements to send to each process
  sendcounts = (int *)malloc(sizeof(int)*size);
  // array describing the displacements where each segment begins
  displs = (int *)malloc(sizeof(int)*size);
  // array describe how many rows sent to every rank
  send_row = (int *)malloc(sizeof(int)*size);
  calc_row = (N-2)/size; //number of rows every rank need to calculta
  extra_row = (N-2)%size; // some ranks need to calculat extra one row elements
  send_element = (calc_row+2+1)*N; // max number of elements some rank receive
  for (i=0; i<size; i++)
  {
    send_row[i] = (i < extra_row) ? calc_row+2+1 : calc_row+2;
    sendcounts[i] = send_row[i] * N;
    displs[i] = begin_segm;
    begin_segm = begin_segm + (sendcounts[i] - 2*N); //each rank displacements start at previous row

  }
  //debug
  if(taskid == MASTER)
  {
    for (i=0; i<size; i++)
    {
      printf("send rank[%d] %d rows, sendcounts[%d] = %d, displs[%d] = %d\n",
    i, send_row[i], i, sendcounts[i], i, displs[i]);
    }
    
    printf("distribute func master rank print initialize matrix A:\n");
    for (i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
      {
        printf("%d\t", A[i*N+j]);
      }
      printf("\n");
    }

  }
  max_element = sendcounts[0];
  rec_buf = (int *)malloc(sizeof(int)*max_element);
  // divide the data among processes as described by sendcounts and displs
  MPI_Scatterv(&A[0], sendcounts, displs, MPI_INT, rec_buf, max_element, MPI_INT, 0, MPI_COMM_WORLD); //use &A[0], not &A(it didnt sent correct value)
  //debug
  printf("rank[%d] receive %d elements:\n", taskid, sendcounts[taskid]);
  for(i=0; i<send_row[taskid]; i++)
  {
    for(j=0; j<N; j++)
    {
      printf("%d\t", rec_buf[i*N+j]);
    }
    printf("\n");
  }
  printf("\n");
  //MPI_Barrier (MPI_COMM_WORLD);
}

void mask_operation (int *A, int N, int*Ap)
{
  int count;
  int k, h;
  printf("mask_operation:rank %d has %d rows, sendcounts[%d] = %d, displs[%d] = %d\n",
  taskid, send_row[taskid],taskid, sendcounts[taskid], taskid, displs[taskid]);
  count = (N-2)*(send_row[taskid]-2); //pixels need be calculated
  result = (int *)malloc(sizeof(int)*count); //allocate space put processed pixels
  int n = 0;
  /*debug*/
  printf("pixels need calculate:\n");
  for(i=1; i<(send_row[taskid]-1); i++)
  {
    for(j=1; j<(N-1); j++)
    {
      h = i*N+j;
      printf("%d\t", rec_buf[h]);
    }
    printf("\n");
  }

  for(i=1; i<(send_row[taskid]-1); i++) // calculate begin second row second column element, without edge
  {
    for(j=1; j<(N-1); j++)
    {
      k = i*N+j; //locate pixel of "e"
      result[n] = (rec_buf[k]*2 + rec_buf[k-N-1]+ rec_buf[k-N]+ rec_buf[k-N+1]+
                   rec_buf[k-1]+rec_buf[k+1]+rec_buf[k+N-1]+rec_buf[k+N]+rec_buf[k+N+1])/10;
      n++;
    }
  }
  
  /*debug*/
  printf("rank %d updated proceesd pixels:\n", taskid);
  for(i=0;i<(send_row[taskid]-2);i++)
  {
    for(j=0;j<(N-2);j++)
    {
      printf("%d\t", result[i*(N-2)+j]);
    }
    printf("\n");
  }
  //MPI_Barrier (MPI_COMM_WORLD);
}

void collect_results (int *Ap, int N)
{
  int rows;
  int sum = 0;
  int numbs = N*N;
  Ap = (int *)malloc(sizeof(int)*numbs); // give processed matrix A space
  send_num = (int *)malloc(sizeof(int)*size);
  /*get elements number that everyrank send to master*/
  for(i=0; i<size; i++)
  {
    if (i == 0 || i == (size-1))
    {
      send_num[i] = sendcounts[i] - N;
    }else{

      send_num[i] = sendcounts[i] - N*2;
    }

    /*allocate displacements for gatherv*/
      displs[i] = sum;
      sum = sum + send_num[i];
  }
  if (taskid == 0 || taskid == (size-1)) //first last rank need send first and last unprocessed row
  {
    rows = send_row[taskid] - 1;
  }else{
    rows = send_row[taskid] - 2;
  }
  send_buf = (int *)malloc(sizeof(int)*send_num[taskid]);
  /*put everyrank edge pixels and processed pixels into send_buf*/
  if (taskid == 0)
  {
    int n=0;
    /*put first edge row to send_buf*/
    for (j=0; j<N;j++)
    {
      send_buf[j] = rec_buf[j];
    }
    /*if column edge, put old value into send_buf, else put processed value*/
    for (i=1; i<rows; i++)
    {
      for(j=0; j<N; j++)
      {
        if(j == 0 || j == (N-1))
        {
          send_buf[i*N+j] = rec_buf[i*N+j];
        }else{
          send_buf[i*N+j] = result[n];
          n++;
        }
      }
    }
  }

  if(taskid == (size-1))
  {
    int n=0;
    int m=0;
    int k= sendcounts[taskid]-N;
    /*put last edge row into send_buf*/
    for (k=sendcounts[taskid]-N; k<sendcounts[taskid];k++)
    {
      m = send_num[taskid]-N+n; //index for last N elements
      send_buf[m] = rec_buf[k];
      n++;
    }
    /*if column edge, put old value into send_buf, else put processed value*/
    n=0;
    for (i=0; i<(rows-1); i++)
    {
      for(j=0; j<N; j++)
      {
        if(j == 0 || j == (N-1))
        {
          send_buf[i*N+j] = rec_buf[(i+1)*N+j];
        }else{
          send_buf[i*N+j] = result[n];
          n++;
        }
      }
    }
  }

  if(taskid > 0 && taskid < (size-1))
  {
    /*if column edge, put old value into send_buf, else put processed value*/
    int n=0;
    for (i=0; i<rows; i++)
    {
      for(j=0; j<N; j++)
      {
        if(j == 0 || j == (N-1))
        {
          send_buf[i*N+j] = rec_buf[(i+1)*N+j];
        }else{
          send_buf[i*N+j] = result[n];
          n++;
        }
      }
    }
  }
  /*debug*/
  printf("rank %d will send send_num: %d, displs: %d, updated sub_matrix is:\n", taskid, send_num[taskid], displs[taskid]);

  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < N; j++)
    {
      printf ("%d\t", send_buf[i*N+j]);
    }
        printf("\n");
  }
  printf("\n");

  MPI_Gatherv(send_buf, send_num[taskid], MPI_INT, Ap, send_num, displs, MPI_INT, 0, MPI_COMM_WORLD);


  //end = clock();
  end = now();
  time_spent = tdiff(begin, end);

  if(taskid == MASTER)
  {
    printf ("total updated Matrix A: \n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf ("%d\t", Ap [i*N+j]);
        }
        printf("\n");
    }
    printf ("\n");
  }
  //MPI_Barrier (MPI_COMM_WORLD);

}
