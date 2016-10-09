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

}
