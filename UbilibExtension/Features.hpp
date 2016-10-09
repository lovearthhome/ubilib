/*
日期:2015.09.30
作者:魏代华,焦帅 余杰
功能：ubilib 特征层 
*/

#ifndef Features_h
#define Features_h
#include<stdlib.h>
#include<math.h>
#include <vector>
#include "FFT.hpp"
using namespace std;
class Features
{
		
public :
    //1求数组最小值
	template<typename T>
	static T miniNum(std::vector<T> data)
	{
		int length = data.size();
		T min = data[0];
		for (int i = 1; i < length; i++)
			min = data[i]<min ? data[i] : min;
		return min;
	}

	//2求数组最大值
	template<typename T>
	static T maxiNum(std::vector<T> data)
	{
		int length = data.size();
		T max = data[0];
		for (int i = 1; i<length; i++)
			max = data[i]<max ? max : data[i];
		return max;
	}

	//3求数组的平均值
	template<typename T>
	static T mean(std::vector<T> data)
	{
		int length = data.size();
		T sum = 0;
		for (int i = 0; i < length; i++)
			sum = sum + data[i];
		return sum / length;
	}

	//4求数组的加权均值
	template<typename T>
	static T weightMean(std::vector<T> data)
	{
		//int length = data.size();
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
	static T standardDeviation(std::vector<T> data)
	{
		int length = data.size();
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
	static T signalVectorMagnitude(std::vector<T> data)
	{
		int length = data.size();
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
	static T median(std::vector<T> data)
	{
		int length = data.size();
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
	static T firstQuartile(std::vector<T> data)
	{
		int length = data.size();
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
	static T thirdQuartile(std::vector<T> data)
	{
		int length = data.size();
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
	static T zeroCrossingRate(std::vector<T> data)
	{
		int length = data.size();
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
	static T meanCrossingsRate(std::vector<T> data)
	{
		int length = data.size();
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
	static T correlationCoefficient(std::vector<T> data1, std::vector<T> data2)
	{
		int length1=data1.size();
		//int length2=data2.size();
		T mean1 = Features::mean(data1);
		T mean2 = Features::mean(data2);
		T covariance = 0.0;
		for (int i = 0; i < length1; i++)
		{
			covariance += (data1[i] - mean1)*(data2[i] - mean2);
		}
		T standarddeviation1 = Features::standardDeviation(data1);
		T standarddeviation2 = Features::standardDeviation(data2);
		return (covariance / length1) / (standarddeviation1*standarddeviation2);
	}

	
	//15求数组的方差
	template<typename T>
	static T variance(std::vector<T> data)
	{
		int length = data.size();
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
	static T spectrumPeakPosition(std::vector<T> data)
	{
		int length = data.size();
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
	static T spectralEnergy(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
		T sum = 0;
		for (int i = 0; i<length; i++){
			sum += pow(data[i], 2);
		}
		return sum / length;
	}

	//18谱熵
	template<typename T>
	static T spectralEntropy(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
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
	static T centroid(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
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
	static T rms(std::vector<T> data){
		int length = data.size();
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
	static T sma(std::vector<T> data, T interval){
		int length = data.size();
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
	static T iqr(std::vector<T> data){
		//int length = data.size();
		return Features::thirdQuartile(data) - Features::firstQuartile(data);
	}

	//23绝对平均值
	template<typename T>
	static T  absMean(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
		T meanD = mean(data);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += fabs(data[i] - meanD);
		}
		return sum / length;
	}

	//24频域偏度
	template<typename T>
	static T  frequencySkew(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
		T meanD = mean(data);
		T dev = standardDeviation(data);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += pow((data[i] - meanD) / dev, 3);
		}
		return sum / length;
	}

	//25频域峰度
	template<typename T>
	static T  frequencyKurt(std::vector<T> data){
		int length = data.size();
		if ( length == 0) return 0;
		T meanD = mean(data);
		T dev = standardDeviation(data);
		T sum = 0;
		int i;
		for (i = 0; i<length; i++){
			sum += pow((data[i] - meanD) / dev, 4) - 3;
		}
		return sum / length;
	}
	
	//26时域能量
	template<typename T>
	static T  tenergy(std::vector<T> data){
		int length = data.size();
		if(length == 0){
			return 0;
		} 
        return spectralEnergy(data);
	}
    
	//27频域标准差
	template<typename T>
	static T  fdev(std::vector<T> data){
		int length = data.size();
		if(length == 0){
			return 0;
		} 
        return standardDeviation(data);
	}
	
	//28频域平均值
	template<typename T>
	static T  fmean(std::vector<T> data){
		int length = data.size();
		if(length == 0){
			return 0;
		} 
        return mean(data);
	}
	
	/*
	 * 因为FFT本身可以对付任意长度的FFT变换，但是为了保险起见，这里还是对傅里叶数据做判断
	 * 本fft不破坏data数组,虽然引用传输
	 * 本fft返回fftdata数组（处理过的一半长度的double数组）
    */
    static std::vector<double>   fft(std::vector<double> &data) {
		std::vector<double> fftdata;
		size_t n = data.size();	//data的实际长度
		unsigned int N = n;//data补齐后的长度	
		if (n < 4 )
		{
			return fftdata;
		}else if ((n & (n - 1)) != 0)  // Is not power of 2,/*如果是2的幂，n一定是100... n-1就是01111....*/
		{
			/*************************************** 
			* unsigned int类型的数二进制中最高位1 
			* 的位置 
			* 例如：0000 0000 0000 0001中返回0 
			*       0000 0000 1111 0000中返回7而不是4 
			*****************************************/ 
			//求出n这个数的二进制的最高的1
			int position=0;  
			unsigned int m=n;  
			while(m)  
			{  
		        m=m>>1;  
		        if(m)  position++;  
		    } 		 
			N = 1<<(position+1);
		}
		//制作real,imag数组，并补齐N
		vector<double> real(n);
		vector<double> imag(n);
		for (vector<double>::iterator it = data.begin(); it != data.end(); ++it)
		{
			real.push_back(*it);
			imag.push_back(0);
		}	
		unsigned int k = N - n;
		while(k--){
			real.push_back(0.0f);
			imag.push_back(0.0f);
		}
		N = real.size();
		if((N & (N - 1)) != 0) 
			throw "fft pre process error";
		else
			FFT::transformRadix2(real, imag);

		//每个fft点，它的模值math.hypot(img,rel)是这个频率上原始信号幅度的（N/2）倍,
		//所以，要把每个fft[i].abs()除以（N/2）才是信号的真正幅度
		unsigned  int fN = N/2;
		for(unsigned int i=1;i<fN+1;i++)
		{
			double amp = std::sqrt(real[i]*real[i] + imag[i]*imag[i]);
			amp = amp * 2.0/(double)N;
			fftdata.push_back(amp);
		}
		//这个返回，fftdata会被c++调用vector的move重新拷贝一个返回，参考
		//http://stackoverflow.com/questions/22655059/why-it-is-ok-to-return-vector-from-function
		return fftdata;
	}	


};
#endif
