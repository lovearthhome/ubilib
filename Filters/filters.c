#include <stdio.h>
#include <stdlib.h>
/**均值滤波**/
double* meanFilter(double *data,int length,int window){
    int i,index;
    double sum;
    double* dataNew=(double*)malloc(length*sizeof(double));
    for(i=0;i<length-window+1;i++){
        sum=0.0;
        for(index=0;index<window;index++){
            if(i==0 || index>window/2){
                dataNew[i+index] = data[i+index];
            }
            sum+=data[i+index];
        }
        dataNew[i+window/2] = sum/window;
    }
    return dataNew;
}

int main()
{
    int i;
    double data[]= {5.2,5.2,5.2,5.2,5.1,25.2,5.2,5.2,5.2,5.2,5.2};
    double* r = meanFilter(data,11,3);
    for(i=0;i<11;i++){
        printf("%d,%.2f\n",i,r[i]);
    }
    printf("\n");
}
