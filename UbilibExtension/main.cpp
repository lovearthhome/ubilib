#include <phpcpp.h>

/**
 *  tell the compiler that the get_module is a pure C function
 */

#include <stdlib.h>
#include "Features.hpp"

using namespace std;
#include <vector>


//因为Php::Value不支持long long,所以要转换一下
Php::Value ubiGetFeatures(Php::Parameters &params)
{
	double samplerate= params[0];
	std::vector<double> accdata;
	Php::Value accdata_php = params[1];
	// assum the value variable holds an array or object, it then
	// is possible to iterate over the values or properties
	for (auto &iter : accdata_php)
    {
		//output key and value
		//下面注释的代码是可以工作的..
	    //Php::out << iter.first << ": " << iter.second<< std::endl;
		double acc = iter.second;
		accdata.push_back(acc);
    }

	//直接使用类名调用成员函数应该使用::，而不是.运算符
	double min = Features::miniNum(accdata);//最小值,y
	double max = Features::maxiNum(accdata);//最大值，y
	double mean = Features::mean(accdata);//均值,y
	double weightMean = Features::weightMean(accdata);//加权均值,n
	double standardDeviation = Features::standardDeviation(accdata);//标准差,y
	double signalVectorMagnitude = Features::signalVectorMagnitude(accdata);//信号向量幅值,n
	double median = Features::median(accdata);//数组的中位数,y
	double firstQuartile = Features::firstQuartile(accdata);//数组1/4分位数,n
	double thirdQuartile = Features::thirdQuartile(accdata);//数组3/4分位数,n
	double zeroCrossingRate = Features::zeroCrossingRate(accdata);//求数组过零率,n
	double meanCrossingsRate = Features::meanCrossingsRate(accdata);//数组过均值率,y
	double correlationCoefficient = Features::correlationCoefficient(accdata,accdata);//相关系数,这个地方只是展示这个函数的使用，本计算没有任何意义在这，将来可以用来计算x轴，y轴之间的相关性,n
	double variance = Features::variance(accdata);//数组的方差,n
	double rms = Features::rms(accdata);//均方根均值,y
	double iqr = Features::iqr(accdata);//四分卫距离,y
	double sma = Features::sma(accdata,1.0f/samplerate);//向量幅值面积,y
	double absMean = Features::absMean(accdata);//绝对平均值,y
	double tenergy = Features::spectralEnergy(accdata);//时域能量,y

	std::vector<double> fftdata = Features::fft(accdata);//FFTdata

	double spectrumPeakPosition = Features::spectrumPeakPosition(fftdata);//数组中的最大值所在的位置-谱峰位置,y
	double spectralEnergy = Features::spectralEnergy(fftdata);//频域能量y
	double spectralEntropy = Features::spectralEntropy(fftdata);//谱熵y
	double centroid = Features::centroid(fftdata);//质心y
	double frequencySkew = Features::frequencySkew(fftdata);//频域偏度y
	double frequencyKurt = Features::frequencyKurt(fftdata);//频域峰度y
	double fdev = Features::standardDeviation(fftdata);//频域标准差y
	double fmean = Features::mean(fftdata);//频域平均值y



	Php::Value retValue;
	retValue["min"] = min;	
	retValue["max"] = max;	
	retValue["mean"] =mean ;	
	retValue["weightMean"] = weightMean;	
	retValue["standardDeviation"] = standardDeviation;	
	retValue["signalVectorMagnitude"] = signalVectorMagnitude;	
	retValue["median"] = median;	
	retValue["firstQuartile"] = firstQuartile;	
	retValue["thirdQuartile"] = thirdQuartile;	
	retValue["zeroCrossingRate"] = zeroCrossingRate;	
	retValue["meanCrossingsRate"] = meanCrossingsRate;	
	retValue["correlationCoefficient"] = correlationCoefficient;	
	retValue["variance"] = variance;	
	retValue["spectrumPeakPosition"] = spectrumPeakPosition;	
	retValue["spectralEnergy"] =spectralEnergy ;	
	retValue["spectralEntropy"] =spectralEntropy ;	
	retValue["centroid"]= centroid;	
	retValue["rms"]= rms;	
	retValue["sma"]=sma ;	
	retValue["iqr"]=iqr ;	
	retValue["absMean"]=absMean ;	
	retValue["frequencySkew"]= frequencySkew;	
	retValue["frequencyKurt"]= frequencyKurt;	
	retValue["fft"]= fftdata;	
	retValue["tenergy"]= tenergy;	
	retValue["fdev"]= fdev;	
	retValue["fmean"]= fmean;	
	return retValue;
}

extern "C" {
    
    /**
     *  Function that is called by PHP right after the PHP process
     *  has started, and that returns an address of an internal PHP
     *  strucure with all the details and features of your extension
     *
     *  @return void*   a pointer to an address that is understood by PHP
     */
    PHPCPP_EXPORT void *get_module() 
    {
        // static(!) Php::Extension object that should stay in memory
        // for the entire duration of the process (that's why it's static)
        static Php::Extension extension("ubilib", "1.0");
		extension.add("ubiGetFeatures", ubiGetFeatures,
			{
				Php::ByVal("samplerate", Php::Type::Float),	//加速度数据的采样频率,	
				Php::ByVal("accdata", Php::Type::Array)//128个float组成的array				
			}); 
        return extension;
    }
}
