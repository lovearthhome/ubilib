#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"
//最小值
double minimum(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double MIN = data[0];
    int i;
    for(i = 1; i < length; i++){
        MIN = data[i]<MIN?data[i]:MIN;
    }
    return MIN;
}
//最大值
double maximum(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double Max = data[0];
    int i;
    for(i = 1; i<length; i++)
        Max = data[i]<Max ? Max : data[i];
    return Max;
}
//方差
double variance(double data[], int length){
    if(data == NULL || length == 0) return 0.0;
    double average = 0, s = 0, sum = 0;
    int i;
    for(i = 0; i<length; i++)
    {
        sum = sum + data[i];
    }
    average = sum / length;
    for(i = 0; i<length; i++)
    {
        s = s + pow(data[i] - average, 2);
    }
    s = s / length;
    return s;
}
//过均值率
double meanCrossingsRate(double data[], int length){
    if(data == NULL || length == 0) return 0.0;
    double Sum = 0;
    double num = 0;
    int i;
    for(i = 0; i < length; i++){
        Sum +=data[i];
    }
    double avg = Sum/length;
    for(i = 0; i < length - 1; i++){
        if ((data[i]-avg) * (data[i + 1]-avg) < 0){
            num++;
        }
    }
    return num / length;
}

//标准差
double standardDeviation(double data[], int length){
    if(data == NULL || length == 0) return 0.0;
    double s = variance(data,length);
    s = sqrt(s);
    return s;
}
//最大值所在的位置=谱峰位置
double spp(double data[],int length){
    if(data == NULL || length == 0) return 0;
    int r = 0;
    double max = data[0];
    int i;
    for(i = 0; i < length; i++) {
        if (data[i] > max ) {
            r = i;
            max = data[i];
        }
    }
    return r;
}
//频域 能量
double energy(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double sum = 0;
    int i;
    for(i=0;i<length;i++){
        sum+=pow(data[i],2);
    }
    return sum/length;
}

//热力学函数  熵
double entropy(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double temp;
    double sum = 0;
    int i;
    for(i=0;i<length;i++){
        temp = data[i];
        if(temp > 0){
            sum+=log(temp)*temp;
        }
    }
    return sum*-1;
}

//平均值
double mean(double data[], int length){
    if(data == NULL || length == 0) return 0.0;

    double Sum = 0;
    int i;
    for(i = 0; i < length; i++)
        Sum = Sum + data[i];
    return Sum / length;
}

// 质心 centroid : SUM(i*f(i)*f(i))/SUM(f(i)*f(i))见Intelligent Sleep Stage Mining Service with Smartphones
double centroid(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double sum1 = 0;
    double sum2 = 0;
    double temp;
    double tempPow;
    int i;
    for(i=0;i<length;i++){
        temp = data[i];
        tempPow = temp*temp;
        sum1+=tempPow;
        sum2+=tempPow*i;
    }
    return sum2/sum1;
}


//向量幅值
double signalVectorMagnitude(double data[], int length){

    double MaxMagnitude = 0;
    int i;
    for(i = 0; i < length; i++){
        if (data[i] < 0)
            MaxMagnitude = (-data[i])<MaxMagnitude ? MaxMagnitude : (-data[i]);
        else
            MaxMagnitude = data[i]<MaxMagnitude ? MaxMagnitude : data[i];
    }
    return MaxMagnitude;
}
//四分位数（Quartile） 1/4
double firstQuartile(double data[], int length){

    if(data == NULL || length == 0) return 0.0;
    double* copydata;
    copydata = (double *)malloc(length*sizeof(double));
    int i,j;
    for(i = 0; i < length; i++){
        copydata[i] = data[i];
    }
    for(i = 0; i<length; i++){
        for (j = i + 1; j<length; j++)
        {
            if (copydata[j]<copydata[i])
            {
                double t;
                t = copydata[i];
                copydata[i] = copydata[j];
                copydata[j] = t;
            }
        }
    }
    double q = 1 + (length - 1) *0.25;
    int b = (int)q;
    double d = q - b;
    double quartile1 = copydata[b-1] + (copydata[b] - copydata[b-1])*d;
    return quartile1;
}
//四分位数（Quartile） 3/4**/
double thirdQuartile(double data[], int length){

    double* copydata;
    copydata = (double*)malloc(length*sizeof(double));
    int i,j;
    for(i = 0; i < length; i++)
    {
        copydata[i] = data[i];
    }
    for(i = 0; i<length; i++)
    {
        for (j = i + 1; j<length; j++)
        {
            if (copydata[j]<copydata[i])
            {
                double t;
                t = copydata[i];
                copydata[i] = copydata[j];
                copydata[j] = t;
            }
        }
    }
    double q = 1 + (length - 1) *0.75;
    int b = (int)q;
    double d = q - b;
    double quartile3 = copydata[b-1] + (copydata[b] - copydata[b-1])*d;
    return quartile3;
}

//过零率
double zeroCrossingRate(double data[], int length){

    double num = 0;
    int i;
    for(i = 0; i < length - 1; i++)
    {
        if (data[i] * data[i + 1]< 0){
            num++;
        }
    }
    return num / length;
}

//中位数（Median）
double median(double data[], int length){

    double* copydata = (double*)malloc(length*sizeof(double));
    int i,j;
    for(i = 0; i < length; i++){
        copydata[i] = data[i];
    }
    double medium;
    for(i = 0; i<length; i++) {
        for (j = i + 1; j<length; j++) {
            if (copydata[j]<copydata[i])
            {
                double t;
                t = copydata[i];
                copydata[i] = copydata[j];
                copydata[j] = t;
            }
        }
    }

    if (length % 2 == 1)
        medium = copydata[length / 2];
    else
        medium = (copydata[length / 2 - 1] + copydata[length / 2]) / 2;
    return medium;
}

/*将某类特征值进行归一化
 * @param lower 归一化的范围最低值
 * @param upper 归一化的范围最高值
 * @param value 要归一化的数据
 * @param min 该类特征的最小值
 * @param max 该类特征的最大值
 * @return
 */
double zeroOneLibSvm(double lower,double upper,double value,double min,double max){
    return lower + (upper-lower) * (value-min)/(max-min);
}
//只能接收长度为128的数组
//若不足128补全为0
double* fftMe(double list[],int length) {
    complex theList[128];
    int i;
    for (i = 0; i < length; i++) {
        complex co;
        co.real=list[i];
        co.imag=0;
        theList[i] = co;
    }
    for (i = length; i < 128; i++) {
        complex co;
        co.imag=0;
        co.real=0;
        theList[i] = co;
    }

    // fft
    fft(length,theList);

    double alpha=1.0/(double)length;
    for (i = 0; i < length; i++) {
        complex co;
        co.real = alpha*theList[i].real;
        co.imag = alpha*theList[i].imag;
        theList[i] = co;
    }
    float fftSeries[length];
    double out[length/2];
    c_abs(theList,fftSeries,length);
    for(i=0;i<length/2;i++){
       out[i] = fftSeries[i+1]*2;
       //printf("%.10f\n",out[i]);
    }
    return out;
}

//均方根平均值
double rms(double list[],int length){
    double rms=0;
    double sum = 0;
    int i;
    for(i=0;i<length;i++){
        sum+=pow(list[i],2);
    }
    rms = sqrt(sum/length);
    return rms;
}

//向量幅值面积 把离散值面积累加起来然后除以总长度。实际是平均每时刻的面积。
double sma(double data[],int length,double interval){
    double sum=0;
    double lot=length * interval;
    int i;
    for(i=0;i<length;i++){
        sum +=data[i]*interval;
    }
    return sum/lot;
}

//四分位距
double iqr(double list[],int length){
    return thirdQuartile(list,length)-firstQuartile(list,length);
}
//绝对平均差
double mad(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double meanD = mean(data,length);
    double sum = 0;
    int i;
    for(i=0;i<length;i++){
        sum+=fabs(data[i]-meanD);
    }
    return sum/length;
}
//时域 能量
double tenergy(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    return energy(data,length);
}
//频域 标准备差
double fdev(double data[],int length){
    return standardDeviation(data,length);
}
//频域平均值
double fmean(double data[],int length){
    return mean(data,length);
}

//频域 偏度
double skew(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double meanD = mean(data,length);
    double dev = standardDeviation(data,length);
    double sum=0;
    int i;
    for(i=0;i<length;i++){
        sum+=pow((data[i]-meanD)/dev,3);
    }
    return sum/length;
}

//频域 峰度
double kurt(double data[],int length){
    if(data == NULL || length == 0) return 0.0;
    double meanD = mean(data,length);
    double dev = standardDeviation(data,length);
    double sum=0;
    int i;
    for(i=0;i<length;i++){
        sum+=pow((data[i]-meanD)/dev,4)-3;
    }
    return sum/length;
}
