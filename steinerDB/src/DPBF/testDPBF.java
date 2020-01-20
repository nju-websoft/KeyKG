package DPBF;

import java.io.IOException;
import java.util.ArrayList;

public class testDPBF {
	GraphDPBF g1;
	String[] filest;
	static String testFile = "E:\\data\\test.txt";
	
	void giveWord() {
		ArrayList<String> aa = new ArrayList<String>();

		aa.add("Владимир");
		aa.add("Gothenburg");	
		aa.add("Hermansverk");
		aa.add("Georgetown");
		g1.givenQueryWord(aa);
	}
	void testDeal() throws IOException {
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "dbpedia201610";
		
		String st = filest[0];		
		g1 = new GraphDPBF();
		g1.Init(st);
				
		DPBF_1 d1 = new DPBF_1();
		
		int iTime= 1;
		double totalnodePruned = 0, weightPruned = 0;
		long PrunedTime = 0;
		long startTime, endTime;
		for (int i = 0; i < iTime; i++) {
			//g1.randQuery(5);			
			giveWord();
			startTime=System.currentTimeMillis();
			totalnodePruned = totalnodePruned + d1.kdeal(g1);
			endTime=System.currentTimeMillis();
			PrunedTime += endTime - startTime;			
			weightPruned += d1.best;						
		}				
		System.out.println("DPBF 平均运行时间： "+PrunedTime/(double)iTime+" ms");
		System.out.println("DPBF 权重： "+ weightPruned/iTime);
		//k1.kdealFile(c1, testFile);
	}
	
	public static void main(String[] args) throws IOException {
		testDPBF t1 = new testDPBF();
		t1.testDeal();	
	}
}
