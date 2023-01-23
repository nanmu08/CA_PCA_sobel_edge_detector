#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char **argv)
{
  /*******begin read image data from txt to array*/
  FILE* readInput;
  readInput = fopen("5000testinput.txt", "r");
  //readInput = fopen("50testinput.txt", "r");
  if (readInput == NULL) {
		printf("File could not opened!");
	}
  int N = 5000;
  //int N = 50;
  int imageArr[N][N];
  int i, j, num;
  for(i = 0; i<N; i++)
  {
    for(j = 0; j<N; j++)
    {
      fscanf(readInput, "%d", &imageArr[i][j]);

    }
  }
  fclose(readInput);
  /*debug*/
  /*printf("read image data array: \n");
  for(i = 0; i<50; i++)
  {
    for(j = 0; j<50; j++)
    {
      printf("%d\t", imageArr[i][j]);
    }
    printf("\n");
  }*/

/********************end reading data to array****************/

/********************begin process data****************/
  int calc_result[N-2][N-2];

  int k, l;
  int Gx = 0, Gy = 0;
  int p1, p2, p3, q1, q2, q3, r1, r2, r3;
  int sumx = 0;
  int sumy = 0;
  int mx[3][3] = {
		{1, 0, -1},
		{2, 0, -2},
		{1, 0, -1}
  };
  int my[3][3] = {
		{1, 2, 1},
		{0, 0, 0},
		{-1, -2, -1}
	};

  for(i=1; i<N-1; i++)
  {
    for(j=1; j<N-1; j++)
    {
      p1 = imageArr[i-1][j-1];
      p2 = imageArr[i-1][j];
      p3 = imageArr[i-1][j+1];
      q1 = imageArr[i][j-1];
      q2 = imageArr[i][j];
      q3 = imageArr[i][j+1];
      r1 = imageArr[i+1][j-1];
      r2 = imageArr[i+1][j];
      r3 = imageArr[i+1][j+1];
      Gx = p1*mx[0][0]+p2*mx[0][1]+p3*mx[0][2]
      +q1*mx[1][0]+q2*mx[1][1]+q3*mx[1][2]
      +r1*mx[2][0]+r2*mx[2][1]+r3*mx[2][2];
      Gy = p1*my[0][0]+p2*my[0][1]+p3*my[0][2]
      +q1*my[1][0]+q2*my[1][1]+q3*my[1][2]
      +r1*my[2][0]+r2*my[2][1]+r3*my[2][2];
      calc_result[i-1][j-1] = sqrt(Gx*Gx) + sqrt(Gy*Gy);
    }
  }
/********************end process data****************/

  /***normalized threshold***/
  int min = 5000, max = 0;

	for(i = 0; i < N-2; i++) {
		for(j = 0; j < N-2; j++) {
			if (calc_result[i][j] < min) {
				min = calc_result[i][j];
			}
			else if (calc_result[i][j] > max) {
				max = calc_result[i][j];
			}
		}
	}

	for(i = 0; i < N-2; i++) {
		for(j = 0; j < N-2; j++) {
			double ratio = (double) (calc_result[i][j] - min) / (max - min);
			calc_result[i][j] = ratio * 255;
		}
	}
  printf("the max value for calc_result is: %d, range is %d\n", max, max-min);


  /****threaholds and denoise****/
  for(i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      if(calc_result[i][j]<51)
      {
        calc_result[i][j] = 0;
      }
    }
  }

  /********debug*******/
  /*printf("after calculated image array is:\n");
  for(i=0; i<N-2; i++)
  {
    for(j=0; j<N-2; j++)
    {
      printf("%d\t", calc_result[i][j]);
    }
    printf("\n");
  }*/

/********************put processed data to outputbuf****************/
  int processedOutput[N][N];
  for (j=0; j<N;j++)
  {
    processedOutput[0][j] = imageArr[0][j];
  }
  for (j=0; j<N;j++)
  {
    processedOutput[N-1][j] = imageArr[N-1][j];
  }
  for (i=1; i<N-1; i++)
    {
      for(j=0; j<N; j++)
      {
        if(j == 0 || j == (N-1))
        {
          processedOutput[i][j] = imageArr[i][j];
        }else{
          processedOutput[i][j] =calc_result[i-1][j-1];
        }
      }
    }
    /*debug*/
  /*printf("whole processed image data array: \n");
  for(i = 0; i<50; i++)
  {
    for(j = 0; j<50; j++)
    {
      printf("%d\t", processedOutput[i][j]);
    }
    printf("\n");
  }*/
  printf("finish print whole processed image data array: \n");
  FILE* writeOutput = fopen("part1Out.txt", "wb");
  //FILE* writeOutput = fopen("50testOut.txt", "wb");
  for(i=0; i<N; i++)
  {
    for(j=0; j<N;j++)
    {
      if(j<N-1)
      {
        fprintf(writeOutput, "%d\t", processedOutput[i][j]);
      }else{
        fprintf(writeOutput, "%d\t\n", processedOutput[i][j]);
      }
    }
  }
  fclose(writeOutput);
}
