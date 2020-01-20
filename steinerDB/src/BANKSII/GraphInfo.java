package BANKSII;

import java.util.*;

import DOGST.Hop;
import database.dbReader;

public class GraphInfo {
	int nodeNum = 0;
	int keywordNum = 0;
	public ArrayList<ArrayList<Hop>> graph;	//read the graph
	public HashMap<String, Integer> HashNodeTable;	//hash Entity Name to Integer
	public ArrayList<String> reHashTable;	//rehash Integer to Entity Name
	public ArrayList<ArrayList<Integer>> invTable;	//keyword inverted table
	public HashMap<String, Integer> HashKeywordTable;	//keyword hash
	public ArrayList<String> reHashKeywordTable;	//rehash Integer to keyword
	public ArrayList<Double> nodeInWeight;
	
	Random r1;
	//read entity hash table
	void readhashTable(String DBName) {
		reHashTable = new ArrayList<String>();
		HashNodeTable = new HashMap<String, Integer>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "subName");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readSubName(a,b)){
			reHashTable.add(a[0]);			
			HashNodeTable.put(a[0], b[0]);
		}
		nodeNum = reHashTable.size();
	}
	
	//read keywrod hash table
	void readKeyMap(String DBName) {
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "keyMap");				
		int []a = new int [1];
		String []b = new String[1];
		while(dh.readKeyMap(a,b)){
			reHashKeywordTable.add(b[0]);
			if (!HashKeywordTable.containsKey(b[0])) {
				HashKeywordTable.put(b[0], a[0]);				
			}
		}
		keywordNum = reHashKeywordTable.size();
	}
	
	//read inverted table
	void readInvTable(String DBName) {
		invTable = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < HashKeywordTable.size(); i++)
			invTable.add(new ArrayList<Integer>());
		
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "invertedTable");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readInvertedTable(a,b)){
			invTable.get(HashKeywordTable.get(a[0])).add(b[0]);
		}
	}
	
	//read graph table
	public void readGraph(String DBName){
		graph = new ArrayList<ArrayList<Hop>>();
		for (int i = 0; i < nodeNum; i++)
			graph.add(new ArrayList<Hop>());
		
		String[] triples = new String[3];
		double[] single = new double[1];
		//input graph
		dbReader dg = new dbReader();
		dg.dbInit(DBName, "graph");
		dbReader dw = new dbReader();
		dw.dbInit(DBName, "edgeWeight");
		while(dg.readGraph(triples) && dw.readEdgeWeight(single)) {
			int source = HashNodeTable.get(triples[0]);
			int target = HashNodeTable.get(triples[2]);
			graph.get(source).add(new Hop(target, single[0]));
			graph.get(target).add(new Hop(source, single[0]));
		}
		
		//calculate incoming reverse edge weight
		nodeInWeight = new ArrayList<Double>();
		for (int i = 0; i < nodeNum; i++) {
			double weight = 0;
			for (Hop it : graph.get(i))
				weight += 1.0/it.dis;
			nodeInWeight.add(weight);
		}
			
		
	}
	
	public void Init(String DBName){	
		readhashTable(DBName);
		readKeyMap(DBName);
		readInvTable(DBName);			
		readGraph(DBName);
		
		r1 = new Random();
		System.out.println("BANKSII init end!");
	}
	
	void output() {
		System.out.println(nodeNum);
	}
	
	public ArrayList<String> keywordList;
	int keyNum;	
	public void randQuery(int kNum){
		keywordList = new ArrayList<String>();		
			
		int kM;
		if (kNum == 0)
			kM = r1.nextInt(6)+2;
		else kM = kNum;
		int gg;
		int []gu = new int[kM];
		for (int i = 0; i < kM; i++) {
			while (true) {
				gg = r1.nextInt(keywordNum);				
				int flag = 0;
				for (int j  = 0; j < i -1 ; j++)	//not generated keyword
					if (gu[j] == gg)
						flag = 1;
				if (flag == 0)
					break;
			}
			keywordList.add(reHashKeywordTable.get(gg));
			gu[i] = gg;
		}
		keyNum = keywordList.size();
	}
	
	public void biggestQuery(int kNum){
		keywordList = new ArrayList<String>();
		PriorityQueue<Hop> nodes = new PriorityQueue<Hop>();
		for (int i = 0; i < invTable.size(); i++)
			nodes.add(new Hop(i, -invTable.get(i).size()));
		
		for (int i = 0; i < kNum; i++) {
			Hop it = nodes.poll();
			System.out.println(it.dis);
			keywordList.add(reHashKeywordTable.get(it.v));
		}
		keyNum = keywordList.size();
	}
	void givenQuery(ArrayList<Integer> a) {
		keywordList = new ArrayList<String>();
		for (int it : a)
			keywordList.add(reHashKeywordTable.get(it));
		keyNum = keywordList.size();
	}
	
	public void givenQueryWord(List<String> a) {
		keywordList = (ArrayList<String>) a;
		keyNum = keywordList.size();
	}
}
