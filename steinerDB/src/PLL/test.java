package PLL;

import java.io.IOException;

import graph.WeightedGraph;

public class test {
	String[] filest;
	test(){
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		filest[2] = "dbpedia201610";
	}
	void testCHL(int loc) throws IOException {				
		String st = filest[loc];		
		WeightedGraph ww = new WeightedGraph();		
		ww.graphDBRead(st);		
		WeightedPLL w2 = new WeightedPLL();
		w2.pllDBDeal(st, ww,"CHL");
	}
	
	void testPLL(int loc) throws IOException {				
		String st = filest[loc];		
		WeightedGraph ww = new WeightedGraph();		
		ww.graphDBRead(st);		
		WeightedPLL w2 = new WeightedPLL();
		w2.pllDBDeal(st, ww,"PLL");	
	}
	public static void main(String[] args) throws IOException {
		test t1 = new test();
		int loc = 0;
		//t1.testCHL(loc);
		t1.testPLL(loc);
	}
}
