package DOGST;

import java.io.IOException;
import BANKSII.*;

public class testDOGST {
	CommonStruct c1;
	String[] filest;
	static String testFile = "E:\\data\\test.txt";
	void testDeal() throws IOException {
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "dbpedia201610";
		
		String st = filest[0];		
		c1 = new CommonStruct();
		c1.Init(st);
				
		DOGST k1 = new DOGST(c1);
		DOGSTSC k2 = new DOGSTSC(c1);
		DOGST k3 = new DOGST(c1);
		
		int iTime= 10;
		double totalnodestar = 0, weightstar = 0;
		double totalnodech = 0, weightchins = 0;
		double totalnodechg = 0, weightchinsg = 0;
		
		long starTime = 0;
		long chinsTime = 0;
		long chinsgTime = 0;	
		long startTime, endTime;
		
		//for bansksII
		GraphInfo g1 = new GraphInfo();
		g1.Init(st);
		BiSearch b1 = new BiSearch();
		double totalnodeBANKS = 0, weightBANKS = 0;
		long banksTime = 0;
		
		for (int i = 0; i < iTime; i++) {
			k1.randQuery(c1, 5);
			/*
			k1.keywordList.clear();
			k1.keywordList.add("Bata");
			k1.keywordList.add("Goyang");
			k1.keywordList.add("Al-Qaseem");
			k1.keywordList.add("Gora");
			*/
			startTime=System.currentTimeMillis();
			totalnodestar = totalnodestar + k1.kdealRandom(c1, k1.keywordList, 0);
			endTime=System.currentTimeMillis();
			starTime += endTime - startTime;
			if (k1.ansList.size()!= 0)
				weightstar += k1.totalweight;
				//weightstar += (k1.totalweight/k1.ansList.size());			
			
			startTime=System.currentTimeMillis();
			totalnodech = totalnodech + k2.kdealRandom(c1, k1.keywordList, 2);
			endTime=System.currentTimeMillis();
			chinsTime += endTime - startTime;
			if (k2.ansList.size()!= 0)
				weightchins += k2.totalweight;
				//weightchins += (k2.totalweight / k2.ansList.size());
			
			startTime=System.currentTimeMillis();
			totalnodechg = totalnodechg + k3.kdealRandom(c1, k1.keywordList, 2);
			endTime=System.currentTimeMillis();
			chinsgTime += endTime - startTime;
			if (k3.ansList.size()!= 0)
				weightchinsg += k3.totalweight;
				//weightchinsg += (k3.totalweight / k3.ansList.size());
			
			g1.givenQueryWord(k1.keywordList);
			startTime=System.currentTimeMillis();
			totalnodeBANKS = totalnodeBANKS + b1.kdealRandom(g1);
			endTime=System.currentTimeMillis();
			banksTime += endTime - startTime;
			if (b1.ansList.size()!= 0)				
				weightBANKS += (b1.totalweight/b1.ansList.size());	
		}				
		System.out.println("star 平均运行时间： "+starTime/(double)iTime+" ms");
		System.out.println("chins 平均运行时间： "+chinsTime/(double)iTime+" ms");
		System.out.println("chinsg 平均运行时间： "+chinsgTime/(double)iTime+" ms");
		System.out.println("banksII 平均运行时间： "+banksTime/(double)iTime+" ms");		
		System.out.println("star 权重： "+ weightstar/iTime);
		System.out.println("chins 权重： "+ weightchins/iTime);		
		System.out.println("chinsg 权重： "+ weightchinsg/iTime);
		System.out.println("banksII 权重： "+ weightBANKS/iTime);
		//k1.kdealFile(c1, testFile);
	}
	public static void main(String[] args) throws IOException {
		testDOGST t1 = new testDOGST();
		t1.testDeal();
		//t1.testDealDisk();	
	}
}
