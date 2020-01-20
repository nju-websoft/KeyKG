package graph;

import java.io.*;
import java.util.*;

import database.dbReader;
//import database.dbWriter;

public class WeightedGraph {
	//static String File = "E:\\data\\exampleU.txt";
	public int nodeNum = 0;
	int keywordNum = 0;
	public ArrayList<ArrayList<Hop>> graph;
	public HashMap<String, Integer> HashNodeTable;	//hash Entity Name to Integer
	public ArrayList<String> reHashTable;	//rehash Integer to Entity Name
	ArrayList<ArrayList<Integer>> invTable;	//keyword inverted table
	HashMap<String, Integer> HashKeywordTable;	//keyword hash
	public ArrayList<String> reHashKeywordTable;	//rehash Integer to keyword
	
	public void graphRead(String testfile) throws IOException{
		File regularFile = new File(testfile+".regular");
		Scanner inputR = new Scanner(regularFile);
		
		File weightFile = new File(testfile+".weight");
		Scanner inputW = new Scanner(weightFile);
		graph = new ArrayList<ArrayList<Hop>>();
		nodeNum = inputR.nextInt();
		for (int i = 0; i < nodeNum; i++)
			graph.add(new ArrayList<Hop>());
		int x, y;
		double w;
		while (inputR.hasNext()) {
			x = inputR.nextInt();
			y = inputR.nextInt();
			w = inputW.nextDouble();
			graph.get(x).add(new Hop(y, w));
			graph.get(y).add(new Hop(x, w));
		}
		inputR.close();
		inputW.close();
	}
	
	//get hash code of st, generate a new hash code if not hashed
	int tryAddVertex(String st) {
		if (!HashNodeTable.containsKey(st)) {
			graph.add(new ArrayList<Hop>());
			HashNodeTable.put(st, nodeNum++); // nodeNum is st's hash number
			reHashTable.add(st);
		}

		return HashNodeTable.get(st);
	}
	
	//read entity hash table
	void readhashTable(String DBName) {
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "subName");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readSubName(a,b)){
			reHashTable.add(a[0]);			
			HashNodeTable.put(a[0], b[0]);
		}
		nodeNum = reHashTable.size();
		for (int i = 0; i < nodeNum; i++)
			graph.add(new ArrayList<Hop>());
	}
		
	public void graphDBRead(String DBName){
		graph = new ArrayList<ArrayList<Hop>>();
		HashNodeTable = new HashMap<String, Integer>();
		reHashTable = new ArrayList<String>();
		String[] triples = new String[3];
		double[] single = new double[1];
		readhashTable(DBName);
		
		//input graph
		dbReader dg = new dbReader();
		dg.dbInit(DBName, "graph");
		dbReader dw = new dbReader();
		dw.dbInit(DBName, "edgeWeight");
		while(dg.readGraph(triples) && dw.readEdgeWeight(single)) {
			//int source = tryAddVertex(triples[0]);
			//int target = tryAddVertex(triples[2]);
			int source = HashNodeTable.get(triples[0]);
			int target = HashNodeTable.get(triples[2]);
			graph.get(source).add(new Hop(target, single[0]));
			graph.get(target).add(new Hop(source, single[0]));
		}		
	}
	
	/*
	//get hash code of st, generate a new hash code if not hashed
	int tryAddKeyword(String st) {
		if (!HashKeywordTable.containsKey(st)) {
			invTable.add(new ArrayList<Integer>());
			HashKeywordTable.put(st, keywordNum++);
			reHashKeywordTable.add(st);
		}

		return HashKeywordTable.get(st);
	}
	
	
	//generate inverted tabel
	public void invertedTableDBGen(String DBName) {
		
		invTable = new ArrayList<ArrayList<Integer>>();
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		dbReader dk = new dbReader();
		dk.dbInit(DBName, "keyword");
		String[] triples = new String[3];	
		while (dk.readKeyword(triples)) 
		if (HashNodeTable.containsKey(triples[0])){
			int source = HashNodeTable.get(triples[0]);
			String[] entity = triples[2].split(" ");
			for (int j = 0; j < entity.length; j++)	//add keyword to all entity's keyList
				invTable.get(tryAddKeyword(entity[j])).add(source); 
		}
		
		//write keyword map
		dbWriter dkm = new dbWriter();
		dkm.dbInit(DBName, "keyMap");
		for (int i = 0; i < keywordNum; i++)
			dkm.insertKeyMap(i, reHashKeywordTable.get(i));
		dkm.dbWriterEnd();
		System.out.println("Keyword Map Generate end!");
		
		dbWriter di = new dbWriter();
		di.dbInit(DBName, "invertedTable");
		for (int i = 0; i < keywordNum; i++)
			for (int j = 0; j < invTable.get(i).size(); j++)
			di.insertInvertedTable(reHashKeywordTable.get(i), invTable.get(i).get(j).intValue());
		di.dbWriterEnd();
		
		System.out.println("Inverted Table Generate end!");
	}
	*/
}
