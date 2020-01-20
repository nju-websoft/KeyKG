package TestCode;

import java.io.IOException;
import java.util.Scanner;

import FileToDB.dbwritepedia;
import PLL.WeightedPLL;
import graph.WeightedGraph;
import graph.WeightAssign;

public class GeneratePLL {
	
	String[] filest;
	
	void testDeal() throws IOException {		
		filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		//filest[2] = "dbpedia201610";
		String cmd1, cmd2, cmd3, cmd4;
		String st = filest[0];
		Scanner sc = new Scanner(System.in);
		System.out.println("Write Graph"+ st +" on database?(Y/N)");		
		cmd1 = sc.nextLine();		
		
		if (cmd1.equals("Y")) {
			System.out.println("Write keyword "+ st +" on database?(Y/N)");			
			cmd2 = sc.nextLine();
			System.out.println("Write invertedTable "+ st +" on database?(Y/N)");			
			cmd3 = sc.nextLine();
			if (st.equals(filest[2])) {
				dbwritepedia pedia = new dbwritepedia();
				pedia.writeGraph("F:\\data\\dbpedia201610\\mappingbased_objects_en.ttl");
				if (cmd2.equals("Y"))
					pedia.writeKeyword("F:\\data\\dbpedia201610\\labels_en.ttl");
				if (cmd3.equals("Y"))
					pedia.invertedTableDBGen(st);				
			}
		}
		System.out.println("write edgeweight " + st + "?(Y/N)");		
		cmd1 = sc.nextLine();
		if (cmd1.equals("Y"))
			WeightAssign.weightDBGen(st);
		System.out.println("Run CHL/PLL on " + st + "?(C/P/N)");		
		cmd4 = sc.nextLine();
		sc.close();
		if (!cmd4.equals("N")) {
			WeightedGraph ww = new WeightedGraph();
			ww.graphDBRead(st);

			WeightedPLL w2 = new WeightedPLL();
			if (cmd4.equals("C"))
				w2.pllDBDeal(st, ww,"CHL");	//this code will change subname and invertedtabel
			if (cmd4.equals("P"))
				w2.pllDBDeal(st, ww, "PLL");
		}
	}
	
	public static void main(String[] args) throws IOException {
		GeneratePLL g1 = new GeneratePLL();
		g1.testDeal();
	}
}
