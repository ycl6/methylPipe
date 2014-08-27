// compile with:
// R CMD SHLIB binning.c
// load in R with:
// dyn.load("binning.so")
#include <R.h>
#include <stdio.h>
//# include <math.h>

void binning(double *score, double *pos, int *poslength, double *binvec, int *binlength, double *binout, int *doAvg) {
  int i,j, newstart;
  double sum, counter;

  newstart = 0;
  for(i=0; i<*binlength; i++) {
    counter=0;
    sum=0;
    for(j=newstart; j<=*poslength; j++) {
     if(pos[j]>binvec[i+1]) {
        newstart = j;
        break;
      }
      if(pos[j]>binvec[i] && pos[j]<=binvec[i+1]) {
        if(score[j]>(-100)) {
          counter++;
          sum=sum+score[j];
        }
      }
    }
    if(counter>0) {
      if(*doAvg>0) {binout[i]=sum/counter;}
      else {binout[i]=sum;}
    }
  }
}


