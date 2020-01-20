package FileToDB;

import java.io.*;
import java.util.*;

import database.dbReader;
import database.dbWriter;
import graph.WeightAssign;


//import com.csvreader.*;

public class dbwritelinkedmdb {
	class Edge implements Comparable<Edge>{
		String e;
		int v;
		
		Edge(String y, int z){
			e = y;
			v =z;
		}
		public int compareTo(Edge e2) {		
			int result = this.v > e2.v ? 1 : (this.v < e2.v ? -1 : 0);			
			return result;
		}
	}
	
	Map<String, Integer> hashNodeTable;
	ArrayList<String> reHashTable;
	ArrayList<Set<Edge>> graph;
	boolean[] visited;
	String database = "linkedmdb";
	
	int tryAdd(String x) {
		if (!hashNodeTable.containsKey(x)) {
			hashNodeTable.put(x, hashNodeTable.size());
			reHashTable.add(x);
			graph.add(new TreeSet<Edge>());
			//System.out.println(hashNodeTable.size());
		}
		return hashNodeTable.get(x);
	}
	
	public void writeGraph(String fileName) throws FileNotFoundException {
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();
		graph = new ArrayList<Set<Edge>>();				
        		
		Scanner input1 = new Scanner(new File(fileName));
		while (input1.hasNext()) {
			String x =  input1.next();
			String y =  input1.next();
			String z =  input1.next();		
			int ux = tryAdd(x);
			int uz = tryAdd(z);
			graph.get(ux).add(new Edge(y, uz));
			graph.get(uz).add(new Edge(y, ux));
		}
		input1.close();
		
		int loc = 0;
		for (int i = 1; i < graph.size(); i++)
			if (graph.get(i).size() > graph.get(loc).size())
				loc = i;
				
		visited = new boolean[graph.size()];
		for (int i = 1; i < graph.size(); i++) visited[i] = false;
		int nodeNum = 0;
		visited[loc] = true;
		ArrayList<Integer> vlist = new ArrayList<Integer>();
		vlist.add(loc);
		while (nodeNum < vlist.size()) {
			int u = vlist.get(nodeNum);
			for (Edge it : graph.get(u))
				if(!visited[it.v]) {
					visited[it.v] = true;
					vlist.add(it.v);
				}
			nodeNum++;
		}
		System.out.println(nodeNum);				
		int cnt = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {				
				cnt++;
		}
		
		int[] edgenum = new int[graph.size()];
		for (int i = 0; i < graph.size(); i++) edgenum[i] = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v])			
						edgenum[i]++;
			}
		
		Map<String, Integer> mapIt = new HashMap<String, Integer>();
		
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						if (!mapIt.containsKey(it.e)) mapIt.put(it.e, 0);
						mapIt.put(it.e, mapIt.get(it.e) + 1);
					}
						
			}
		
		
		dbWriter dkm = new dbWriter();
		dkm.dbInit(database, "graph"); 
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						dkm.insertGraph(reHashTable.get(i), it.e, reHashTable.get(it.v));
						edgenum[i]++;
					}
			}
		dkm.dbWriterEnd();
		System.out.println("Graph edges write end!");
		
		dbWriter dw = new dbWriter();
		dw.dbInit(database, "edgeWeight");
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v])		
						dw.insertEdgeWeight((Math.log(mapIt.get(it.e)+1)/Math.log(2.0)));					
			}
		dw.dbWriterEnd();			
		System.out.println("Edge weights write end!");
		
		
		cnt = 0;
		dbWriter d1= new dbWriter();
		d1.dbInit(database, "subName");
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {				
				d1.insertSubName(reHashTable.get(i), cnt);
				cnt++;
			}
		d1.dbWriterEnd();
		System.out.println("Sub Name Old write end!");
	}
	
	public void writeKeyword(String fileName) throws FileNotFoundException {
		
		dbWriter dw = new dbWriter();
		dw.dbInit(database, "keyword");
		Scanner input1 = new Scanner(new File(fileName));		
		while (input1.hasNext()) {
			String[] line = input1.nextLine().split("	");
			if(line.length != 2) continue;
			String x =  line[0];		
			String z =  line[1];			
			z = z.replace("\"", "");			
			
			if (hashNodeTable.containsKey(x))
				if (visited[hashNodeTable.get(x)])
					dw.insertKeyword(x, "0", z);
			
		}
		input1.close();
		dw.dbWriterEnd();		
	}
	
	int tryAddKeyword(String x) {
		if (!HashKeywordTable.containsKey(x)) {
			invTable.add(new ArrayList<Integer>());
			HashKeywordTable.put(x, HashKeywordTable.size());
			reHashKeywordTable.add(x);
		}
		return HashKeywordTable.get(x);
			
	}
	Map<String, Integer>HashKeywordTable;
	ArrayList<String>reHashKeywordTable;
	ArrayList<ArrayList<Integer>>invTable;
	public void invertedTableDBGen() {
		
		hashNodeTable = new HashMap<String, Integer>();
		invTable = new ArrayList<ArrayList<Integer>>();
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		
		dbReader dsb = new dbReader();
		dsb.dbInit(database, "subName");			
		String []a = new String[1];
		int []b = new int [1];
		while(dsb.readSubName(a,b))
			hashNodeTable.put(a[0], b[0]);		
		
		dbReader dk = new dbReader();
		dk.dbInit(database, "keyword");
		String[] triples = new String[3];
		while (dk.readKeyword(triples)) 
		if (hashNodeTable.containsKey(triples[0])){
			int source = hashNodeTable.get(triples[0]);
			triples[2]= triples[2].toLowerCase();
			triples[2] = triples[2].replace("(", " ");
			triples[2] = triples[2].replace(")", " ");
			triples[2] = triples[2].replaceAll(" +", " ");
			String[] entity = triples[2].split(" ");
			for (int j = 0; j < entity.length; j++)	//add keyword to all entity's keyList
				invTable.get(tryAddKeyword(entity[j])).add(source); 
		}
		
		//write keyword map
		dbWriter dkm = new dbWriter();
		dkm.dbInit(database, "keyMap");
		for (int i = 0; i < reHashKeywordTable.size(); i++)
			dkm.insertKeyMap(i, reHashKeywordTable.get(i));
		dkm.dbWriterEnd();
		System.out.println("Keyword Map Generate end!");
		
		dbWriter di = new dbWriter();
		di.dbInit(database, "invertedTable");
		for (int i = 0; i < reHashKeywordTable.size(); i++)
			for (int j : invTable.get(i))
			di.insertInvertedTable(reHashKeywordTable.get(i), j);
		di.dbWriterEnd();
		
		System.out.println("Inverted Table Generate end!");
	}

	
	void dbdeal() throws IOException  {
		String []fileName = new String[3]; 
				
		fileName[0] = "F:\\data\\linkedmdb\\linkedmdb.txt";
		//writeGraph(fileName[0]);
		
		fileName[1] = "F:\\data\\linkedmdb\\keyword.txt";
		//writeKeyword(fileName[1]);

		WeightAssign.weightDBGen(database);
		//invertedTableDBGen();
	}
	
	public static void main(String[] args) throws IOException {
		dbwritelinkedmdb link = new dbwritelinkedmdb();
		link.dbdeal();		
	}
}
