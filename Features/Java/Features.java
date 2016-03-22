package com.ydcun.ivita.util;
import java.util.Arrays;

public class Features{

	/**求数组最小值 double**/
	public static double minimum(double data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		double MIN = data[0];
		for (int i = 1; i < length; i++){
			MIN = data[i]<MIN?data[i]:MIN;
		}
		return MIN;
	}
	/**求数组最小值 float**/
	public static double minimum(float data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		float MIN = data[0];
		for (int i = 1; i < length; i++){
			MIN = data[i]<MIN?data[i]:MIN;
		}
		return MIN;
	}
	
	/**求数组最小值 int**/
	public static double minimum(int data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		double MIN = data[0];
		for (int i = 1; i < length; i++){
			MIN = data[i]<MIN?data[i]:MIN;
		}
		return MIN;
	}

	/**求数组最大值 double**/
	public static double maximum(double data[]){
		if(data == null || data.length == 0) return 0.0;
		
		int length = data.length;
		double Max = data[0];
		for (int i = 1; i<length; i++)
			Max = data[i]<Max ? Max : data[i];
		return Max;
	}
	
	/**求数组最大值 float**/
	public static double maximum(float data[]){
		if(data == null || data.length == 0) return 0.0;		
		int length = data.length;
		float Max = data[0];
		for (int i = 1; i<length; i++)
			Max = data[i]<Max ? Max : data[i];
		return Max;
	}
	
	/**求数组最大值 int **/
	public static double maximum(int data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		int Max = data[0];
		for (int i = 1; i<length; i++)
			Max = data[i]<Max ? Max : data[i];
		return Max;
	}

	/**求数组的方差**/
	public static double variance(double data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		double average = 0, s = 0, sum = 0;
		for (int i = 0; i<length; i++)
		{
			sum = sum + data[i];
		}
		average = sum / length;
		for (int i = 0; i<length; i++)
		{
			s = s + Math.pow(data[i] - average, 2);
		}
		s = s / length;
		return s;
	}

	/**求数组的方差 float**/
	public static double variance(float data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		double average = 0, s = 0, sum = 0;
		for (int i = 0; i<length; i++)
		{
			sum = sum + data[i];
		}
		average = sum / length;
		for (int i = 0; i<length; i++)
		{
			s = s + Math.pow(data[i] - average, 2);
		}
		s = s / length;
		return s;
	}

	/**求数组的方差 int **/
	public static double variance(int data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		double average = 0, s = 0, sum = 0;
		for (int i = 0; i<length; i++)
		{
			sum = sum + data[i];
		}
		average = sum / length;
		for (int i = 0; i<length; i++)
		{
			s = s + Math.pow(data[i] - average, 2);
		}
		s = s / length; 
		return s;
	}

	/**求double数组过均值率**/
	public static double meanCrossingsRate(double data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		double Sum = 0;
		double num = 0;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
			Sum +=copydata[i];
		}
		double avg = Sum/length;
		for (int i = 0; i < length; i++)
		{
			copydata[i] = copydata[i] - avg;
		}
		for (int i = 0; i < length - 1; i++)
		{
			if (copydata[i] * copydata[i + 1]< 0){
				num++;
			}
		}
		return num / length;
	}
	
	/**求float数组过均值率**/
	public static double meanCrossingsRate(float data[]){
		if(data == null || data.length == 0) return 0.0;	
		int length = data.length;
		float Sum = 0;
		float num = 0;
		float[] copydata = new float[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
			Sum +=copydata[i];
		}
		float avg = Sum/length;
		for (int i = 0; i < length; i++)
		{
			copydata[i] = copydata[i] - avg;
		}
		for (int i = 0; i < length - 1; i++)
		{
			if (copydata[i] * copydata[i + 1]< 0){
				num++;
			}
		}
		return num / length;
	}


	
	/**求double数组的标准差**/
	public static double standardDeviation(double data[]){
		if(data == null || data.length == 0) return 0.0;	
		double s = variance(data);
		s = Math.sqrt(s); 
		return s;
	}

	/**求float数组的标准差 float**/
	public static double standardDeviation(float data[]){
		if(data == null || data.length == 0) return 0.0;	
		double s = variance(data);
		s = Math.sqrt(s); 
		return s;
	}


	/**
	 * 求一个double数组中的最大值所在的位置=谱峰位置
	 * @param data FFT array
	 * @return
	 */
	public static double spp(double[] data){
		if(data == null || data.length == 0) return 0;	
		int r = 0;
		double max = data[0];
		for (int i = 0; i < data.length; i++) {
			if (data[i] > max ) {
				r = i;
				max = data[i];
			}
		}
		return r;
	}
	
	
	/**
	 * 频域 能量
	 * FIXME:如果data数组有误，应该返回什么？
	 * @param data FFT array
	 */
	public static double energy(double[] data){
		if(data == null || data.length == 0) return 0.0;	
		double sum = 0;
		for(int i=0;i<data.length;i++){
			sum+=Math.pow(data[i],2);
		}
		return sum/data.length;
	}
	
	/**
	 * 热力学函数  熵
	 * java.lang.Math.log(double a) 返回自然对数（以e为底）的一个double值。特殊情况：如果参数是NaN或小于零，那么结果是NaN.如果参数是正无穷大，那么结果是正无穷大。如果参数是正零或负零，那么结果是负无穷大。
	 * @param data  有loge=ln运算不允许有和零值
	 * @return
	 * @throws InfoException 
	 */
	public static double entropy(double[] data){
		if(data == null || data.length == 0) return 0.0;
		double temp;
		double sum = 0;
		for(int i=0;i<data.length;i++){
			temp = data[i];
			if(temp > 0){
				sum+=Math.log(temp)*temp;
			}
		}
		return sum*-1;
	}

	/**求数组的平均值 double**/
	public static double mean(double data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		double Sum = 0;
		for (int i = 0; i < length; i++)
			Sum = Sum + data[i];
		return Sum / length;
	}
	
	/**求数组的平均值 float**/
	public static double mean(float data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		float Sum = 0;
		for (int i = 0; i < length; i++)
			Sum = Sum + data[i];
		return Sum / length;
	}

	/**
	 * 质心 centroid : SUM(i*f(i)*f(i))/SUM(f(i)*f(i))见Intelligent Sleep Stage Mining Service with Smartphones 
	 * @param data
	 * @return
	 */
	public static double centroid(double[] data){
		if(data == null || data.length == 0) return 0.0;
		double sum1 = 0;
		double sum2 = 0;
		double temp;
		double tempPow;
		for(int i=0;i<data.length;i++){
			temp = data[i];
			tempPow = temp*temp;
			sum1+=tempPow;
			sum2+=tempPow*i;
		}
		return sum2/sum1;
	}
	
	
	/**求信号向量幅值 （Signal vector magnitude）**/
	public static double signalVectorMagnitude(double data[]){
		int length = data.length;
		double MaxMagnitude = 0;
		for (int i = 0; i < length; i++){
			if (data[i] < 0)
				MaxMagnitude = (-data[i])<MaxMagnitude ? MaxMagnitude : (-data[i]);
			else
				MaxMagnitude = data[i]<MaxMagnitude ? MaxMagnitude : data[i];
		}
		return MaxMagnitude;
	}
	
	
	/**求信号向量幅值 （Signal vector magnitude）**/
	public static double signalVectorMagnitude(float data[]){
		int length = data.length;
		double MaxMagnitude = 0;
		for (int i = 0; i < length; i++){
			if (data[i] < 0)
				MaxMagnitude = (-data[i])<MaxMagnitude ? MaxMagnitude : (-data[i]);
			else
				MaxMagnitude = data[i]<MaxMagnitude ? MaxMagnitude : data[i];
		}
		return MaxMagnitude;
	}
	/**求信号向量幅值 （Signal vector magnitude）**/
	public static double signalVectorMagnitude(int data[]){
		int length = data.length;
		double MaxMagnitude = 0;
		for (int i = 0; i < length; i++){
			if (data[i] < 0)
				MaxMagnitude = (-data[i])<MaxMagnitude ? MaxMagnitude : (-data[i]);
			else
				MaxMagnitude = data[i]<MaxMagnitude ? MaxMagnitude : data[i];
		}
		return MaxMagnitude;
	}
	
	
	/**
	 * 
	 * 四分位数（Quartile） 1/4
	 * 
	 * FIXME: 这里使用了一个低效的排序算法,需要改进。
	 * 
	 * **/
	public static double firstQuartile(double data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++){
			copydata[i] = data[i];
		}
		for (int i = 0; i<length; i++){
			for (int j = i + 1; j<length; j++)
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
	
	/**
	 * 
	 * 四分位数（Quartile） 1/4
	 * 
	 * FIXME: 这里使用了一个低效的排序算法,需要改进。
	 * **/
	public static double firstQuartile(float data[]){
		if(data == null || data.length == 0) return 0.0;
		int length = data.length;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++){
			copydata[i] = data[i];
		}
		for (int i = 0; i<length; i++){
			for (int j = i + 1; j<length; j++)
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

	/**四分位数（Quartile） 3/4**/
	public static double thirdQuartile(double data[]){
		int length = data.length;
		double[] copydata = new double[length];
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
	
	/**四分位数（Quartile） 3/4**/
	public static double thirdQuartile(float data[]){
		int length = data.length;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		for (int i = 0; i<length; i++)
		for (int j = i + 1; j<length; j++)
		{
			if (copydata[j]<copydata[i])
			{
					double t;
					t = copydata[i];
					copydata[i] = copydata[j];
					copydata[j] = t;
			}
		}

		double q = 1 + (length - 1) *0.75;
		int b = (int)q;
		double d = q - b;
		double quartile3 = copydata[b-1] + (copydata[b] - copydata[b-1])*d;
		return quartile3;
	}
	/**四分位数（Quartile） 3/4**/
	public static double thirdQuartile(int data[]){
		int length = data.length;
		double[] copydata = new double[length];
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
	/**求数组过零率**/
	public static double zeroCrossingRate(double data[]){
		int length = data.length;
		double num = 0;
		for (int i = 0; i < length - 1; i++)
		{
			if (data[i] * data[i + 1]< 0){
				num++;
			}
		}
		return num / length;
	}
	/**求数组过零率**/
	public static double zeroCrossingRate(float data[]){
		int length = data.length;
		double num = 0;
		for (int i = 0; i < length - 1; i++)
		{
			if (data[i] * data[i + 1]< 0){
				num++;
			}
		}
		return num / length;
	}

	/**相关系数*/
	public static double correlation(double data1[], double data2[]){
		int length1 = data1.length;
		int length2 = data2.length;
		double mean1 = mean(data1);
		double mean2 = mean(data2);
		double covariance = 0.0;
		for (int i = 0; i < length1; i++)
		{
			covariance += (data1[i] - mean1)*(data2[i] - mean2);
		}
		double standarddeviation1 = standardDeviation(data1);
		double standarddeviation2 = standardDeviation(data2);
		return (covariance / length1) / (standarddeviation1*standarddeviation2);
	}
	/**相关系数*/
	public static double correlation(float data1[], float data2[]){
		int length1 = data1.length;
		int length2 = data2.length;
		double mean1 = mean(data1);
		double mean2 = mean(data2);
		double covariance = 0.0;
		for (int i = 0; i < length1; i++)
		{
			covariance += (data1[i] - mean1)*(data2[i] - mean2);
		}
		double standarddeviation1 = standardDeviation(data1);
		double standarddeviation2 = standardDeviation(data2);
		return (covariance / length1) / (standarddeviation1*standarddeviation2);
	}


	/**中位数（Median）**/  
	public static double median(double data[]){
		int length = data.length;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		double medium;
		for (int i = 0; i<length; i++)
		{
			for (int j = i + 1; j<length; j++)
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

		if (length % 2 == 1)
			medium = copydata[length / 2];
		else
			medium = (copydata[length / 2 - 1] + copydata[length / 2]) / 2;
		return medium;
	}
	/**中位数（Median）**/  
	public static double median(float data[]){
		int length = data.length;
		double[] copydata = new double[length];
		for (int i = 0; i < length; i++)
		{
			copydata[i] = data[i];
		}
		double medium;
		for (int i = 0; i<length; i++)
		{
			for (int j = i + 1; j<length; j++)
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
		if (length % 2 == 1)
			medium = copydata[length / 2];
		else
			medium = (copydata[length / 2 - 1] + copydata[length / 2]) / 2;
		return medium;
	}

	/**特征规范化最值法
	 * 同一列的每个特征除以该列特征中的最大值将特征范围变为（-1,1）
	 * */
	public static double[][] minmax(double[][] data){
		int r = data.length;
		int c = data[0].length;
		double[][] temp =new double[r][c];
		for (int i = 0; i<r; i++)
		{
			temp[i] = new double[c];
		}
		//求每一列的最大值和最小值
		double[][] middata = new double[r][c];
		middata[0] = new double[c];
		middata[1] = new double[c];
		for (int i = 0; i<c; i++)
		{
			double midmin = data[0][i];
			double midmax = data[0][i];
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
			double a = middata[1][i] - middata[0][i];
			for (int j = 0; j<r; j++)
			{
				double b = data[j][i] - middata[0][i];
				temp[j][i] = b / a;
			}
		}
		return temp;
	}

	/**特征规范化0-1法
	 * 同一列的每个特征减去该列特征的平均值，然后除以该列特征的标准方差
	 * */
	public static double[][] zeroOne(double[][] data){
		int r = data.length;
		int c = data[0].length;
		double[][] temp;/*声明一个二维数组，用来存储处理过的数据*/
		temp = new double[r][c];
		for (int i = 0; i<r; i++){
			temp[i] = new double[c];
		}
		/*求每一列的平均值和标准方差*/
		double[][] middata;
		middata = new double[r][c];
		middata[0] = new double[c];
		middata[1] = new double[c];
		for (int i = 0; i<c; i++){
			double[] mid = new double[r];
			for (int j = 0; j<r; j++)
			{
				mid[j] = data[j][i];
			}
			double midmean = mean(mid);
			double midstandarddeviation = standardDeviation(mid);
			middata[0][i] = midmean;
			middata[1][i] = midstandarddeviation;
		}
		//0-1法处理数据，存储在数组temp中返回
		for (int i = 0; i<c; i++){
			for (int j = 0; j<r; j++)
			{
				double b = data[j][i] - middata[0][i];
				temp[j][i] = b / middata[1][i];
			}
		}
		return temp;
	}
	/**
	 * 将某类特征值进行归一化
	 * @param lower 归一化的范围最低值
	 * @param upper 归一化的范围最高值
	 * @param value 要归一化的数据
	 * @param min 该类特征的最小值
	 * @param max 该类特征的最大值
	 * @return
	 */
	public static double zeroOneLibSvm(double lower,double upper,double value,double min,double max){
		return lower + (upper-lower) * (value-min)/(max-min);
	}
    /*只能接收长度为128的数组 
    若不足128补全为0
    */
    public static double[] fft(double[] list) {
		int len=list.length;
		Complex[] theList = new Complex[128];
		for (int i = 0; i < 128; i++) {
			theList[i] = new Complex(list[i], 0);
		}
		for (int i = len; i < 128; i++) {
			theList[i] = new Complex(0, 0);
		}
		
		// fft
		Complex[] Y = FFT.fft(theList);

		double alpha=1.0/(double)len;
		for (int i = 0; i < Y.length; i++) {
			Y[i] = Y[i].times(alpha);
		}	
		double[] fftSeries = new double[64];
		for (int i = 1, j = 0; i < 64 + 1; i++, j++) {
			fftSeries[j] = 2 * Y[i].abs();
		}
	
		return fftSeries;
	}	
    
    /**
     * 均方根平均值
     * @param list
     * @return
     */
    public static double rms(double[] list){
    	double rms=0;
    	double sum = 0;
    	for(int i=0;i<list.length;i++){
    		sum+=Math.pow(list[i],2);
    	}
    	rms = Math.sqrt(sum/list.length);
    	return rms;
    }
    
    /**
     * 向量幅值面积
     * 把离散值面积累加起来然后除以总长度。实际是平均每时刻的面积。
     * 最后*1/T
     * @param list
     * @return
     */
    public static double sma(double[] data,double interval){
    	double sum=0;
    	double lot=data.length * interval;
    	for(int i=0;i<data.length;i++){
    		sum +=data[i]*interval;
    	}
    	return sum/lot;
    }
    
    /**
     * 四分卫距
     * @param list
     * @return
     */
    public static double iqr(double[] list){
    	return Features.thirdQuartile(list)-Features.firstQuartile(list);
    }
    /**
     * 绝对平均差
     * @param list
     * @return
     */
    public static double mad(double[] data){
		if(data == null || data.length == 0) return 0.0;
    	double mean = Features.mean(data);
    	double sum = 0;
    	for(int i=0;i<data.length;i++){
    		sum+=Math.abs(data[i]-mean);
    	}
    	return sum/data.length;
    }
    /**
     * 时域 能量
     * @param list
     * @return
     */
    public static double tenergy(double[] data){
		if(data == null || data.length == 0) return 0.0;
    	return Features.energy(data);
    }
    /**
     * 频域 标准备差
     * @param list
     * @return
     */
    public static double fdev(double[] data){
    	return Features.standardDeviation(data);
    }
    /**
     * 频域平均值
     * @param list
     * @return
     */
    public static double fmean(double[] data){
    	return Features.mean(data);
    }
   
    /**
     * 频域 偏度
     * @param list
     * @return
     */
    public static double skew(double[] data){
		if(data == null || data.length == 0) return 0.0;
    	double mean = Features.mean(data);
    	double dev = Features.standardDeviation(data);
    	double sum=0;
    	for(int i=0;i<data.length;i++){
    		sum+=Math.pow((data[i]-mean)/dev,3);
    	}
    	return sum/data.length;
    }
    
    /**
     * 频域 峰度
     * @param list
     * @return
     */
    public static double kurt(double[] data){
		if(data == null || data.length == 0) return 0.0;
    	double mean = Features.mean(data);
    	double dev = Features.standardDeviation(data);
    	double sum=0;
    	for(int i=0;i<data.length;i++){
    		sum+=Math.pow((data[i]-mean)/dev,4)-3;
    	}
    	return sum/data.length;
    }
    
}
