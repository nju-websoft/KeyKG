package graph;

import java.io.*;
import java.util.*;

import database.dbReader;
import database.dbWriter;

public class WeightAssign {

	//static String testFile = "E:\\data\\exampleU.txt";
	static String testFile = "E:\\data\\undirectedGraph(dbpedia201510).txt";
	public void weightGen(String testfile) throws IOException{

		Random ran=new Random();
		File regularFile = new File(testfile+".regular");
		Scanner input = new Scanner(regularFile);
		File weightFile = new File(testfile + ".weight");
		if (!weightFile.exists())
			weightFile.createNewFile();
		PrintWriter weightdata = new PrintWriter(weightFile);
		input.nextInt();
		while (input.hasNextInt()) {
			input.nextInt();
			input.nextInt();
			weightdata.println(ran.nextFloat()*1000);
			//weightdata.println(1);
		}
		input.close();
		weightdata.close();
	}
	
	public static void weightDBGen(String DBName) {
		Scanner s1 = new Scanner(System.in);
		System.out.println("Dangerous operation! Overwrite edgeweight on\n" + DBName +"\n?(YY/N)");
		String choise = s1.nextLine();
		while (!choise.equals("N") && !choise.equals("YY")) {
			System.out.println("Input 'YY' or 'N'");
			choise = s1.nextLine();			
		}
		if (choise.equals("N")) {
			s1.close();
			return;
		}
		s1.close();
		
		dbReader dg = new dbReader();
		dg.dbInit(DBName, "graph");
		dbWriter dw = new dbWriter();
		
		dw.dbInit(DBName, "edgeWeight");
		Random ran=new Random();
		String[] triples = new String[3];		
		while(dg.readGraph(triples)) {
			dw.insertEdgeWeight(ran.nextFloat()*99+1);			
		}		
		dw.dbWriterEnd();
		System.out.println("Edge Weight Generate end!");
	}
	
	public void DBtoFile() throws FileNotFoundException {
		String []filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		filest[2] = "dbpedia201610";		
		
		String DBName = filest[1];
		String outputFile = "E:\\data\\" + DBName + "DB.txt";
		PrintWriter pw = new PrintWriter(new File(outputFile));
		String[] triples = new String[3];
		double[] single = new double[1];
		
		Map<String, Integer> HashNodeTable = new TreeMap<>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "subName");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readSubName(a,b))			
			HashNodeTable.put(a[0], b[0]);		

		dbReader dg = new dbReader();
		dg.dbInit(DBName, "graph");
		dbReader dw = new dbReader();
		dw.dbInit(DBName, "edgeWeight");
		while(dg.readGraph(triples) && dw.readEdgeWeight(single)) {
			int source = HashNodeTable.get(triples[0]);
			int target = HashNodeTable.get(triples[2]);
			//int len = (int) (single[0]*100);
			double len = single[0];
			pw.println(source + " " + target + " " + len);
			//pw.println(source + " " + target);
		}	
		
		pw.close();
	}
	public static void main(String[] args) throws IOException {
		WeightAssign W1 = new WeightAssign();
		W1.DBtoFile();
		//W1.weightGen(testFile);
	}
}
