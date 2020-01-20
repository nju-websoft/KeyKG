package prunedDP;

import java.io.IOException;
import java.util.ArrayList;

public class testPrunedDP {
	GraphDP g1;
	String[] filest;
	static String testFile = "E:\\data\\test.txt";
	
	void giveWord() {
		ArrayList<String> aa = new ArrayList<String>();
		aa.add("carey");
		aa.add("great");
		aa.add("gatsby");		
		aa.add("shame");
		aa.add("education");
		aa.add("never");	
		aa.add("let");	
		aa.add("go");	
		aa.add("greatest");	
		g1.givenQueryWord(aa);
	}
	void testDeal() throws IOException {
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		filest[2] = "dbpedia201610";
		
		String st = filest[1];		
		g1 = new GraphDP();
		g1.Init(st);
				
		//PrunedNaive p1 = new PrunedNaive();
		PrunedBiHeap p1 = new PrunedBiHeap();
		
		int iTime= 1;
		double totalnodePruned = 0, weightPruned = 0;
		long PrunedTime = 0;
		long startTime, endTime;
		for (int i = 0; i < iTime; i++) {
			//g1.randQuery(5);			
			giveWord();
			startTime=System.currentTimeMillis();
			totalnodePruned = totalnodePruned + p1.kdeal(g1);
			endTime=System.currentTimeMillis();
			PrunedTime += endTime - startTime;			
			weightPruned += p1.best;						
		}				
		System.out.println("Pruned++ 平均运行时间： "+PrunedTime/(double)iTime+" ms");
		System.out.println("Pruned++ 权重： "+ weightPruned/iTime);
		//k1.kdealFile(c1, testFile);
	}
	
	public static void main(String[] args) throws IOException {
		testPrunedDP t1 = new testPrunedDP();
		t1.testDeal();	
	}
}
