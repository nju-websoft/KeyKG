package DPBF;

import java.io.IOException;
import java.util.*;

import database.dbReader;

public class GraphDPBF {
	int nodeNum = 0;
	int keywordNum = 0;
	public ArrayList<ArrayList<EdgeDPBF>> graph;	//read the graph
	public HashMap<String, Integer> HashNodeTable;	//hash Entity Name to Integer
	public ArrayList<String> reHashTable;	//rehash Integer to Entity Name
	public ArrayList<ArrayList<Integer>> invTable;	//keyword inverted table
	public HashMap<String, Integer> HashKeywordTable;	//keyword hash
	public ArrayList<String> reHashKeywordTable;	//rehash Integer to keyword
	//public ArrayList<Double> nodeInWeight;
	
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
		graph = new ArrayList<ArrayList<EdgeDPBF>>();
		for (int i = 0; i < nodeNum; i++)
			graph.add(new ArrayList<EdgeDPBF>());
		
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
			graph.get(source).add(new EdgeDPBF(source, target, single[0]));
			graph.get(target).add(new EdgeDPBF(target, source, single[0]));
		}
		
	}
	
	public void Init(String DBName){	
		readhashTable(DBName);
		readKeyMap(DBName);
		readInvTable(DBName);			
		readGraph(DBName);
		
		System.out.println("DPBF init end!");
	}
	
	
	void output() {
		System.out.println(nodeNum);
	}
	
	ArrayList<String> keywordList;
	int keyNum;	
	void randQuery(int kNum){
		keywordList = new ArrayList<String>();		
		Random r1 = new Random();	
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
	
	void givenQuery(ArrayList<Integer> a) {
		keywordList = new ArrayList<String>();
		for (int it : a)
			keywordList.add(reHashKeywordTable.get(it));
		keyNum = keywordList.size();
	}
	
	public void givenQueryWord(ArrayList<String> a) {
		keywordList = a;
		keyNum = keywordList.size();
	}
	
	
	void testDeal() {
		Init("mondial");
		Scanner sc = new Scanner(System.in);
		double ans = 0;
		while(sc.hasNext()) {
			int u = sc.nextInt();
			int v = sc.nextInt();
			boolean noedge = true;
			for (EdgeDPBF it : graph.get(u))
				if (it.v == v) {
					ans = it.weight;
					noedge = false;
				}
			if (noedge) {
				System.out.println("No such edges!");
				break;
			}
			System.out.println(ans);
		}
		sc.close();
	}
	public static void main(String[] args) throws IOException {
		GraphDPBF g1 = new GraphDPBF();
		g1.testDeal();	
	}
}
