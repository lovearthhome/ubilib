/*
日期:2015.09.30
作者:魏代华,焦帅
功能：ubilib 特征层 
*/

#ifndef Features_h
#define Features_h
#include<stdlib.h>
#include<math.h>

class Features
{
		
public :
    //1求数组最小值
	template<typename T>
	static T miniNum(T data[], int length)
	{
		T min = data[0];
		for (int i = 1; i < length; i++)
			min = data[i]<min ? data[i] : min;
		return min;
	}

	//2求数组最大值
	template<typename T>
	static T maxiNum(T data[], int length)
	{
		T max = data[0];
		for (int i = 1; i<length; i++)
			max = data[i]<max ? max : data[i];
		return max;
	}

	//3求数组的平均值
	template<typename T>
	static T mean(T data[], int length)
	{
		T sum = 0;
		for (int i = 0; i < length; i++)
			sum = sum + data[i];
		return sum / length;
	}

	//4求数组的加权均值
	template<typename T>
	static T weightMean(T data[])
	{
		double sum = 0;
		T result;
		double w[5] = { 0.15, 0.15, 0.2, 0.2, 0.3 };
		for (int i = 0; i < 5; i++)
			sum = sum + data[i + 5] * w[i];
		result = T(sum);
		return result;
	}

	//5求数组的标准差
	template<typename T>
	static T standardDeviation(T data[], int length)
	{
		T average = 0, s = 0, sum = 0;
		for (int i = 0; i<length; i++)
		{
			sum = sum + data[i];
		}
		average = sum / length;
		for (int i = 0; i<length; i++)
		{
			s = s + pow(data[i] - average, 2);// 偏离平均数的距离和 
		}
		s = s / length;//方差 
		s = sqrt(s); //标准差 
		return s;
	}

	//6求信号向量幅值
	template<typename T>
	static T signalVectorMagnitude(T data[], int length)
	{
		T maxMagnitude = 0;
		for (int i = 0; i < length; i++){
			if (data[i] < 0)
				maxMagnitude = (-data[i])<maxMagnitude ? maxMagnitude : (-data[i]);
			else
				maxMagnitude = data[i]<maxMagnitude ? maxMagnitude : data[i];
		}
		return maxMagnitude;
	}

	//7求数组的中位数
	template<typename T>
	static T median(T data[], int length)
	{
		T *copydata = new T[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		T medium;
		for (int i = 0; i<length; i++)
		{
			for (int j = i + 1; j<length; j++)
			{
				if (copydata[j]<copydata[i])
				{
					T t;
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

	//8求数组1/4分位数
	template<typename T>
	static T firstQuartile(T data[], int length)
	{
		T *copydata = new T[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		for (int i = 0; i<length; i++)
		{
			for (int j = i + 1; j<length; j++)
			{
				if (copydata[j]<copydata[i])
				{
					T t;
					t = copydata[i];
					copydata[i] = copydata[j];
					copydata[j] = t;
				}
			}
		}
		T q = 1 + (length - 1) *0.25;
		int b = int(q);
		T d = q - b;
		T quartile1 = copydata[b-1] + (copydata[b] - copydata[b-1])*d;
		return quartile1;
	}

	//9求数组3/4分位数
	template<typename T>
	static T thirdQuartile(T data[], int length)
	{
		T *copydata = new T[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		for (int i = 0; i<length; i++)
		{
			for (int j = i + 1; j<length; j++)
			{
				if (copydata[j]<copydata[i])
				{
					T t;
					t = copydata[i];
					copydata[i] = copydata[j];
					copydata[j] = t;
				}
			}
		}
		T q = 1 + (length - 1) *0.75;
		int b = int(q);
		T d = q - b;
		T quartile3 = copydata[b-1] + (copydata[b] - copydata[b-1])*d;
		return quartile3;
	}

	//10求数组过零率
	template<typename T>
	static T zeroCrossingRate(T data[], int length)
	{
		T num = 0;
		for (int i = 0; i < length - 1; i++)
		{
			if (data[i] * data[i + 1]< 0){
				num++;
			}
		}
		return num / length;
	}

	//11求数组过均值率
	template<typename T>
	static T meanCrossingsRate(T data[], int length)
	{
		T sum = 0;
		int num = 0;
		T *copydata = new T[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
			sum = sum + data[i];
		}
		T avg = sum / length;
		for (int i = 0; i < length; i++)
		{
			copydata[i] = copydata[i] - avg;
		}
		for (int i = 0; i < length - 1; i++)
		{
			if (copydata[i] * copydata[i+1]< 0){
				num++;
			}
		}
		return (T)num / length;
	}

	//12求相关系数
	template<typename T>
	static T correlationCoefficient(T data1[], int length1, T data2[], int length2)
	{
		T mean1 = Features::mean(data1, length1);
		T mean2 = Features::mean(data2, length2);
		T covariance = 0.0;
		for (int i = 0; i < length1; i++)
		{
			covariance += (data1[i] - mean1)*(data2[i] - mean2);
		}
		T standarddeviation1 = Features::standardDeviation(data1, length1);
		T standarddeviation2 = Features::standardDeviation(data2, length2);
		return (covariance / length1) / (standarddeviation1*standarddeviation2);
	}

	//13特征规范化最值法
	template<typename T>
	static T ** minMax(T data[][5], int r, int c)
	{
		T **temp;//声明一个二维数组，用来存储处理过的数据
		temp = new T*[c];
		for (int i = 0; i<r; i++)
		{
			temp[i] = new T[c];
		}
		//求每一列的最大值和最小值
		T **middata;
		middata = new T*[c];
		middata[0] = new T[c];
		middata[1] = new T[c];
		for (int i = 0; i<c; i++)
		{
			T midmin = data[0][i];
			T midmax = data[0][i];
			for (int j = 0; j<r; j++)
			{
				if (midmin>data[j][i])
				{
					midmin = data[j][i];
				}
				if (midmax<data[j][i])
				{
					midmax = data[j][i];
				}
			}
			middata[0][i] = midmin;
			middata[1][i] = midmax;
		}
		//最值法处理数据，存储在数组temp中返回
		for (int i = 0; i<c; i++)
		{
			T a = middata[1][i] - middata[0][i];
			for (int j = 0; j<r; j++)
			{
				T b = data[j][i] - middata[0][i];
				temp[j][i] = b / a;
			}
		}
		return temp;
	}

	//14特征规范化0-1法
	template<typename T>
	static T ** zeroOne(T data[][5], int r, int c)
	{
		T **temp;//声明一个二维数组，用来存储处理过的数据
		temp = new T*[c];
		for (int i = 0; i<r; i++)
		{
			temp[i] = new T[c];
		}
		//求每一列的平均值和标准方差
		T **middata;
		middata = new T*[c];
		middata[0] = new T[c];
		middata[1] = new T[c];
		for (int i = 0; i<c; i++)
		{
			T *mid = new T[r];
			for (int j = 0; j<r; j++)
			{
				mid[j] = data[j][i];
			}
			T midmean = Features::mean(mid, r);
			T midstandarddeviation = Features::standardDeviation(mid, r);
			middata[0][i] = midmean;
			middata[1][i] = midstandarddeviation;
		}
		//0-1法处理数据，存储在数组temp中返回
		for (int i = 0; i<c; i++)
		{
			for (int j = 0; j<r; j++)
			{
				T b = data[j][i] - middata[0][i];
				temp[j][i] = b / middata[1][i];
			}
		}
		return temp;
	}

	//15求数组的方差
	template<typename T>
	static T variance(T data[], int length)
	{
		if (length == 0) return 0;
		T average = 0, s = 0, sum = 0;
		for (int i = 0; i<length; i++)
		{
			sum = sum + data[i];
		}
		average = sum / length;
		for (int i = 0; i<length; i++)
		{
			s = s + pow(data[i] - average, 2);// 偏离平均数的距离和 
		}
		s = s / length;//方差 
		return s;
	}

	//16数组中的最大值所在的位置-谱峰位置
	template<typename T>
	static T spectrumPeakPosition(T data[], int length)
	{
		if (length == 0)return 0;
		T max = data[0];
		int location = 0;
		for (int i = 1; i < length; i++){
			if (data[i] > max) {
				location = i;
				max = data[i];
			}
		}
		return location;
	}

	//17频域能量
	template<typename T>
	static T spectralEnergy(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T sum = 0;
		for (int i = 0; i<length; i++){
			sum += pow(data[i], 2);
		}
		return sum / length;
	}

	//18谱熵
	template<typename T>
	static T spectralEntropy(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T temp;
		T sum = 0;
		for (int i = 0; i<length; i++){
			temp = data[i];
			if (temp > 0){
				sum += log(temp)*temp;
			}
		}
		return sum*-1;
	}

	//19质心
	template<typename T>
	static T centroid(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T sum1 = 0;
		T sum2 = 0;
		T temp;
		T tempPow;
		for (int i = 0; i<length; i++){
			temp = data[i];
			tempPow = temp*temp;
			sum1 += tempPow;
			sum2 += tempPow*i;
		}
		return sum2 / sum1;
	}

	//20均方根均值
	template<typename T>
	static T rms(T data[], int length){
		T rms = 0;
		T sum = 0;
		for (int i = 0; i<length; i++){
			sum += pow(data[i], 2);
		}
		rms = sqrt(sum / length);
		return rms;
	}

	//21向量幅值面积 把离散值面积累加起来然后除以总长度。实际是平均每时刻的面积。
	template<typename T>
	static T sma(T data[], int length, T interval){
		T sum = 0;
		T lot = length * interval;
		int i;
		for (i = 0; i<length; i++){
			sum += data[i] * interval;
		}
		return sum / lot;
	}

	//22四分卫距离
	template<typename T>
	static T iqr(T data[], int length){
		return Features::thirdQuartile(data, length) - Features::firstQuartile(data, length);
	}

	//23绝对平均值
	template<typename T>
	static T  absMean(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T meanD = mean(data, length);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += fabs(data[i] - meanD);
		}
		return sum / length;
	}

	//24频域偏度
	template<typename T>
	static T  frequencySkew(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T meanD = mean(data, length);
		T dev = standardDeviation(data, length);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += pow((data[i] - meanD) / dev, 3);
		}
		return sum / length;
	}

	//25频域峰度
	template<typename T>
	static T  frequencyKurt(T data[], int length){
		if (data == NULL || length == 0) return 0;
		T meanD = mean(data, length);
		T dev = standardDeviation(data, length);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += pow((data[i] - meanD) / dev, 4) - 3;
		}
		return sum / length;
	}

};
#endif
