package FileToDB;

import java.io.*;
import java.util.*;

import org.apache.jena.rdf.model.*;
import org.apache.jena.util.FileManager;

import database.dbReader;
import database.dbWriter;


//import com.csvreader.*;

public class dbwritepedia {
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
	String database = "dbpedia201610";
	Random r1;
	int tryAdd(String x) {
		if (!hashNodeTable.containsKey(x)) {
			hashNodeTable.put(x, hashNodeTable.size());
			reHashTable.add(x);
			graph.add(new TreeSet<Edge>());
			System.out.println(hashNodeTable.size());
		}
		return hashNodeTable.get(x);
	}
	
	public void writeGraph(String fileName) {
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();
		graph = new ArrayList<Set<Edge>>();		
		Model model = ModelFactory.createDefaultModel();
        InputStream in = FileManager.get().open( fileName );
        model.read(in, null, "TTL");            
        StmtIterator iter = model.listStatements();
        
		String x, y, z;
		
		while (iter.hasNext()) {
			org.apache.jena.rdf.model.Statement stmt = iter.nextStatement(); // get next statement
			//RDFNode object = stmt.getObject();
			x = stmt.getSubject().toString(); // get the subject
			y = stmt.getPredicate().toString(); // get the predicate
			z = stmt.getObject().toString();			
			int ux = tryAdd(x);
			int uz = tryAdd(z);
			graph.get(ux).add(new Edge(y, uz));
			graph.get(uz).add(new Edge(y, ux));
		}
				
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
		while(nodeNum < vlist.size()) {
			int u = vlist.get(nodeNum);
			for (Edge it : graph.get(u))
				if(!visited[it.v]) {
					visited[it.v] = true;
					vlist.add(it.v);
				}
			nodeNum++;
		}
		System.out.println(nodeNum);
		
		
		int[] edgenum = new int[graph.size()];
		
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
		
		r1 = new Random();
		dbWriter dw = new dbWriter();
		dw.dbInit(database, "edgeWeight");
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						double nn = edgenum[i];
						if (edgenum[it.v] > nn) nn = edgenum[it.v];
						nn = Math.log(1 + nn)/Math.log(2.0);
						//dw.insertEdgeWeight(nn);
						dw.insertEdgeWeight(r1.nextDouble()*100);
					}
			}
		dw.dbWriterEnd();
		System.out.println("Edge weights write end!");
		
		int cnt = 0;
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
	
	public void writeKeyword(String fileName) {
		Model model = ModelFactory.createDefaultModel();
        InputStream in = FileManager.get().open( fileName );
        model.read(in, null, "TTL");            
        StmtIterator iter = model.listStatements();
        
		String x, y, z;
		dbWriter dw = new dbWriter();
		dw.dbInit(database, "keyword");
		while (iter.hasNext()) {
			org.apache.jena.rdf.model.Statement stmt = iter.nextStatement(); // get next statement
			//RDFNode object = stmt.getObject();
			x = stmt.getSubject().toString(); // get the subject
			y = stmt.getPredicate().toString(); // get the predicate
			z = stmt.getObject().toString();
			z = z.replace("@en","");
			//if (object instanceof Resource) 
			if (hashNodeTable.containsKey(x))
				if (visited[hashNodeTable.get(x)])
					dw.insertKeyword(x, y, z);
		}
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
	public void invertedTableDBGen(String DBName) {
		
		hashNodeTable = new HashMap<String, Integer>();
		invTable = new ArrayList<ArrayList<Integer>>();
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		
		dbReader dsb = new dbReader();
		dsb.dbInit(DBName, "subName");			
		String []a = new String[1];
		int []b = new int [1];
		while(dsb.readSubName(a,b))
			hashNodeTable.put(a[0], b[0]);		
		
		dbReader dk = new dbReader();
		dk.dbInit(DBName, "keyword");
		String[] triples = new String[3];	
		while (dk.readKeyword(triples)) 
		if (hashNodeTable.containsKey(triples[0])){
			int source = hashNodeTable.get(triples[0]);
			triples[2]= triples[2].toLowerCase();
			String[] entity = triples[2].split(" ");
			for (int j = 0; j < entity.length; j++)	//add keyword to all entity's keyList
				invTable.get(tryAddKeyword(entity[j])).add(source); 
		}
		
		//write keyword map
		dbWriter dkm = new dbWriter();
		dkm.dbInit(DBName, "keyMap");
		for (int i = 0; i < reHashKeywordTable.size(); i++)
			dkm.insertKeyMap(i, reHashKeywordTable.get(i));
		dkm.dbWriterEnd();
		System.out.println("Keyword Map Generate end!");
		
		dbWriter di = new dbWriter();
		di.dbInit(DBName, "invertedTable");
		for (int i = 0; i < reHashKeywordTable.size(); i++)
			for (int j : invTable.get(i))
			di.insertInvertedTable(reHashKeywordTable.get(i), j);
		di.dbWriterEnd();
		
		System.out.println("Inverted Table Generate end!");
	}
	
	String outfile = "E:\\data\\dbpedia201610.txt";;
	public void writeToOut(String fileName) throws IOException {
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();
		graph = new ArrayList<Set<Edge>>();		
		Model model = ModelFactory.createDefaultModel();
        InputStream in = FileManager.get().open( fileName );
        model.read(in, null, "TTL");            
        StmtIterator iter = model.listStatements();
        
		String x, y, z;
		
		while (iter.hasNext()) {
			org.apache.jena.rdf.model.Statement stmt = iter.nextStatement(); // get next statement
			//RDFNode object = stmt.getObject();
			x = stmt.getSubject().toString(); // get the subject
			y = stmt.getPredicate().toString(); // get the predicate
			RDFNode ob = stmt.getObject();			
			z = ob.toString();						
			if(!(ob instanceof Resource)) {
				System.out.println(z);
			}
			
			int ux = tryAdd(x);
			int uz = tryAdd(z);
			graph.get(ux).add(new Edge(y, uz));
			graph.get(uz).add(new Edge(y, ux));
		}
		//if (true) return;
				
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
		
		
		
		Map<String, Integer> mapIt = new HashMap<String, Integer>();
		
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						if (!mapIt.containsKey(it.e)) mapIt.put(it.e, 0);
						mapIt.put(it.e, mapIt.get(it.e) + 1);
					}
						
			}
		
			
		
		int[] edgenum = new int[graph.size()];
		for (int i = 0; i < graph.size(); i++) edgenum[i] = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v])			
						edgenum[i]++;
			}
		
		int name[] = new int[graph.size()];
		int cnt = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				name[i] = cnt;
				cnt++;
		}
		
		PrintWriter fopinfo;
		File file = new File(outfile);
		if (!file.exists()) file.createNewFile();
		fopinfo = new PrintWriter(file);
		fopinfo.println(cnt);
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						double nn = edgenum[i];
						if (edgenum[it.v] > nn) nn = edgenum[it.v];
						nn = Math.log(1 + nn)/Math.log(2.0);
						fopinfo.println(name[i] + " " + name[it.v]);
						//fopinfo.println(name[i] + " " + name[it.v] + " " + nn);
						//fopinfo.println(name[i] + " " + name[it.v] + " " + (int)(Math.log(mapIt.get(it.e))/Math.log(2.0)));
					}
			}
		fopinfo.close();
		System.out.println("Edge weights write end!");
	}
	
	void dbdeal() throws IOException  {
		String []fileName = new String[3]; 
				
		fileName[0] = "F:\\data\\dbpedia201610\\mappingbased_objects_en.ttl";
		writeToOut(fileName[0]);

	}
	
	public static void main(String[] args) throws IOException {
		dbwritepedia link = new dbwritepedia();
		link.dbdeal();		
	}
}
