package com.ydcun.ivita.util;
import java.util.Arrays;

public class Filters{

	/*MedianFilter author MingyuYi*/
	//支持偶数滤波，当滤波窗口是偶数的时候，L选大，R选小。如窗口大小为4时，L=2，R=1
	//当filter窗口是偶数时候要以中间的两个数的平均值作为中位数。参考如下链接
	//http://jingyan.baidu.com/article/425e69e69f07fbbe15fc161f.html
	public static double[] medianFilter(double[] list,int len){
		int L=0,R=0;        
		if(len%2==0){
			L=len/2;
			R=len/2-1;
		}else{
			L=R=len/2;
		}
		double[] data=new double[list.length];
		double[] filter=new double[len];
		for(int i=0;i<list.length;i++){
			for(int j=0;j<len;j++){
				int p=i+j-L;
				filter[j]=((p>=0)&&(p<list.length))?list[p]:0;
			}
			Arrays.sort(filter);
			double mid = 0.0f;
			if(len%2 == 0)
			{
				mid = (filter[len/2-1]+filter[len/2])/2.0;
			}else{
				
				mid = filter[len/2];
			}
			data[i]=mid;
		}
		return data;
	}
  
	/*BandpassFilter author MingyuYi*/
	//只让某个频段的频率通过。从参考中看到，ifft后只需要实数部分就可以了
	/*
	 jAVA参考：http://introcs.cs.princeton.edu/java/97data/FFT.java.html
   /***************************************************************************
    *  Test client and sample execution
    *
    *  % java FFT 4
    *  x
    *  -------------------
    *  -0.03480425839330703
    *  0.07910192950176387
    *  0.7233322451735928
    *  0.1659819820667019
    *
    *  y = fft(x)
    *  -------------------
    *  0.9336118983487516
    *  -0.7581365035668999 + 0.08688005256493803i
    *  0.44344407521182005
    *  -0.7581365035668999 - 0.08688005256493803i
    *
    *  z = ifft(y)
    *  -------------------
    *  -0.03480425839330703
    *  0.07910192950176387 + 2.6599344570851287E-18i
    *  0.7233322451735928
    *  0.1659819820667019 - 2.6599344570851287E-18i
    *
    *
	
	*/
	public static double[] bandpassFilter(double[] list,int minp,int maxp){
		int len=list.length;
		Complex[] theList = new Complex[128];
		for (int i = 0; i < len; i++) {
			theList[i] = new Complex(list[i], 0);
		}
		for (int i = len; i < 128; i++) {
			theList[i] = new Complex(0, 0);
		}
		
		// fft
		Complex[] Y = FFT.fft(theList);
		//注意，fft频谱是对称的，要清理左右
		//必须保证0分量不变，否则会丢失整个信号的直流分量
		//举例，0-1-2-3-4-5-6-7
		//0是直流分量，长度N=8，minp=1,maxp=2,那么 6 = 8 - maxp  7 = N - minp
		int N = Y.length;
		for (int i = 1; i < Y.length; i++) {
			if(i>=minp && i<=maxp || i>=N - maxp && i<=N - minp)
			{
				//do nothing
			}else{
				Y[i] = new Complex(0, 0);
			}	
		}
		//
		Complex[] X = FFT.ifft(Y);
		for (int i = 0; i < X.length; i++) {
			list[i] = X[i].re();
		}
		return list;
	}

	public static double[] sppsFilter(double[] list,int K/*频谱排序的大小*/){
		int len=list.length;
		Complex[] theList = new Complex[128];
		for (int i = 0; i < len; i++) {
			theList[i] = new Complex(list[i], 0);
		}
		for (int i = len; i < 128; i++) {
			theList[i] = new Complex(0, 0);
		}
		
		// fft
		Complex[] Y = FFT.fft(theList);
		int N = Y.length;
		double[] filter = new double[N/2];
		//注意，fft频谱是对称的，要清理左右
		//必须保证0分量不变，否则会丢失整个信号的直流分量
		//举例，0-1-2-3-4-5-6-7
		//0是直流分量，长度N=8，minp=1,maxp=2,那么 6 = 8 - maxp  7 = N - minp
		filter[0] = 0.0f;
		for (int i = 0; i < N/2; i++) 
		{
			filter[i] = Y[i].abs(); 
		}
		/*
		 * 理应，还需要对fft的结果 
		 * 实部虚部除以 N
		 * 绝对值乘以  2
		 * 但我们仅仅是为了获取排序关系，所以不需要这么做了
		 * 注意：arrays.sort是从小到大排序: 
		 * */
		Arrays.sort(filter);	

		if(K >= N/2 ){K = N/2 -1;}
		if(K <   0  ){K =  0;    }

		double threshold = filter[N/2 - K];//因为是从小到大排序，所以需要.

		for (int i = 1; i < N; i++) {
			if(Y[i].abs() < threshold)
			{
				Y[i] = new Complex(0, 0);
			}
		}

		Complex[] X = FFT.ifft(Y);
		for (int i = 0; i < X.length; i++) {
			list[i] = X[i].re();
		}
		return list;
	}
	
	/**
	 * JKalman 例子
	 * 下载JKalman  https://sourceforge.net/projects/jkalman/
     * Main method
     * @param args
     */
    public static void main(String[] args) {
    	double[] data = {11.171030601657474,9.904405622953956,8.939142405092579,8.939142405092579,8.30643354866536,8.30643354866536,10.943998785555719,10.943998785555719,10.943998785555719,12.203577935742242,12.49899734844664,12.49899734844664,10.679458629195736,10.679458629195736,9.092961195006483,9.092961195006483,8.262315551852094,8.262315551852094,9.262648859344937,10.802380473828755,10.802380473828755,13.693880305129733,13.693880305129733,13.693880305129733,13.693880305129733,15.091629358598665,15.091629358598665,10.0141731202998,10.0141731202998,8.93955655917932,8.93955655917932,7.495405127651238,7.495405127651238,8.922286289735487,8.922286289735487,11.551180977219278,11.551180977219278,11.551180977219278,14.980502711848812,13.409765014006279,13.409765014006279,13.409765014006279,10.864727772749164,10.864727772749164,8.093524940102652,8.093524940102652,7.92466499691148,7.92466499691148,9.168958285849003,9.168958285849003,11.915573275880528,11.915573275880528,11.915573275880528,11.915573275880528,15.157452467174487,15.157452467174487,13.189046790082122,13.189046790082122,10.973476131748933,10.973476131748933,7.733405493124161,7.733405493124161,7.625072269632799,7.625072269632799,9.678337800227403,9.678337800227403,12.61587194820492,12.61587194820492,12.61587194820492,14.770776807983827,10.997423513218704,10.997423513218704,10.997423513218704,9.487827065294027,8.783601917887394,8.783601917887394,8.197584608991567,8.197584608991567,10.895354212785344,10.895354212785344,10.895354212785344,10.895354212785344,12.34316051963155,12.34316051963155,13.536513178688004,13.536513178688004,11.10770643432613,11.10770643432613,9.357416647182264,9.357416647182264,7.827420873193128,7.827420873193128,8.568715807071369,8.568715807071369,10.57140403954836,10.57140403954836,12.679919262674773,12.679919262674773,12.679919262674773,14.557761749673329,10.251802324037417,10.251802324037417,10.251802324037417,9.006324094517641,8.28957899715677,8.28957899715677,8.122806324848655,8.122806324848655,8.122806324848655,8.122806324848655,11.449475954615572,11.449475954615572,12.687444668151718,12.687444668151718,14.678078050588887,14.678078050588887,9.930578768314914,9.930578768314914,9.18400106355575,9.18400106355575,8.024783221501783,8.024783221501783,9.036688347378274,9.036688347378274,10.801377593442718,10.801377593442718,13.733648532082356,13.733648532082356};
        try {
            JKalman kalman = new JKalman(1, 1);
            double[][] A= new double[][]{{1}};
            double[][] H= new double[][]{{1}};
            double[][] Q= new double[][]{{1}};
            double[][] R= new double[][]{{6}};
            kalman.setTransition_matrix(new Matrix(A));
            kalman.setMeasurement_matrix(new Matrix(H));
            kalman.setProcess_noise_cov(new Matrix(Q));
            kalman.setMeasurement_noise_cov(new Matrix(R));
            kalman.setError_cov_post(kalman.getError_cov_post().identity());
            
            //开始位置
            Matrix statePost = new Matrix(1,1);
            statePost.set(0, 0, data[0]);
            kalman.setState_post(statePost);
            
            Matrix measurementZ = new Matrix(1,1);
            Matrix predictX = null;
            Matrix currectionX = null;
            System.out.print(data[0]+" "+statePost);
            for(int i=1;i<128;i++){
            	measurementZ.set(0, 0,data[i]);
            	predictX = kalman.Predict();
            	currectionX = kalman.Correct(measurementZ);
            	System.out.print(data[i] +" "+currectionX);
            }
            
            
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
    }

}
