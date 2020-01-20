package BANKSII;

import java.io.*;
import java.util.*;


public class testBanks {
	GraphInfo g1;
	String[] filest;
	static String testFile = "E:\\data\\test.txt";
	
	void giveWord() {
		ArrayList<String> aa = new ArrayList<String>();				
		aa.add("eiffel");
		g1.givenQueryWord(aa);
	}
	void testDeal() throws IOException {
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "dbpedia201610";
		
		String st = filest[1];		
		g1 = new GraphInfo();
		g1.Init(st);
				
		BiSearch b1 = new BiSearch();
		
		int iTime= 1;
		double totalnodeBANKS = 0, weightBANKS = 0;
		long banksTime = 0;
		long startTime, endTime;
		for (int i = 0; i < iTime; i++) {
			//g1.randQuery(5);			
			giveWord();
			startTime=System.currentTimeMillis();
			totalnodeBANKS = totalnodeBANKS + b1.kdealRandom(g1);
			endTime=System.currentTimeMillis();
			banksTime += endTime - startTime;
			if (b1.ansList.size()!= 0)
				weightBANKS += b1.totalweight/b1.ansList.size();						
		}				
		System.out.println("banksII 平均运行时间： "+banksTime/(double)iTime+" ms");
		System.out.println("banksII 权重： "+ weightBANKS/iTime);
		//k1.kdealFile(c1, testFile);
	}
	
	public static void main(String[] args) throws IOException {
		testBanks t1 = new testBanks();
		t1.testDeal();	
	}
}
