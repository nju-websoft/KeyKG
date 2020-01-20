package FileToDB;

import java.io.*;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import org.apache.jena.rdf.model.*;
import org.apache.jena.rdf.model.Statement;
import org.apache.jena.util.FileManager;

import FileToDB.dbwritepedia.Edge;
import graph.WeightAssign;
//import com.csvreader.*;

public class dbwritemondial {

	static final String JDBC_DRIVER = "com.mysql.jdbc.Driver";
	static final String DB_URL = "jdbc:mysql://localhost:3306/?user=root&characterEncoding=utf8";
	static final String USER = "root";
	static final String PASS = "shiyuxuan";
	
	void writeAll(String fileName, String database, String[] table) {
		Connection conn = null;
		java.sql.Statement stmt = null;
		PreparedStatement pstm0, pstm1, pstm2;
		try {
			Class.forName("com.mysql.jdbc.Driver");
			System.out.println("Connecting to database...");
			conn = DriverManager.getConnection(DB_URL, USER, PASS);
			stmt = conn.createStatement();
			String sql0, sql1 , sql2;
			sql0 = "CREATE DATABASE IF NOT EXISTS " + database;
			stmt.executeUpdate(sql0);
			sql0 = "USE " + database;
			stmt.executeUpdate(sql0);
			for (int i = 0; i < table.length; i++) {
				sql0 = "DROP TABLE IF EXISTS " + table[i];
				stmt.executeUpdate(sql0);
				sql0 = "CREATE TABLE IF NOT EXISTS " + table[i] + "(subject VARCHAR(1023)," + "predicate VARCHAR(1023),"
						+ " object VARCHAR(1023))";
				stmt.executeUpdate(sql0);
			}
			
			conn.setAutoCommit(false);
			sql0 = "INSERT INTO " + table[0] + "(subject, predicate, object) VALUES (?,?,?)";
			pstm0 = conn.prepareStatement(sql0);
			
			sql1 = "INSERT INTO " + table[1] + "(subject, predicate, object) VALUES (?,?,?)";
			pstm1 = conn.prepareStatement(sql1);
			
			sql2 = "INSERT INTO " + table[2] + "(subject, predicate, object) VALUES (?,?,?)";
			pstm2 = conn.prepareStatement(sql2);

			Model model = ModelFactory.createDefaultModel();
			InputStream in = FileManager.get().open(fileName);
			model.read(in, null, "RDF/XML-ABBREV");
			
			String x, y, z;

			ResIterator subjects = model.listSubjects();
			int i = 0;
			while (subjects.hasNext()) {
				Resource subject = subjects.next();
				// get all triples
				StmtIterator properties = subject.listProperties();
				while (properties.hasNext()) {
					Statement stmtS = properties.nextStatement();					
					RDFNode object = stmtS.getObject();
					x = subject.toString();
					y = stmtS.getPredicate().toString(); // get the predicate
					z = object.toString();
					if (object instanceof Resource) 
					{						
						pstm0.setString(1, x);
						pstm0.setString(2, y);
						pstm0.setString(3, z);
						pstm0.addBatch();
						i++;						
					}		
					if (y.indexOf("name") != -1) {
						pstm1.setString(1, x);
						pstm1.setString(2, y);
						pstm1.setString(3, z);
						pstm1.addBatch();
						i++;
					}
					if (y.indexOf("type") != -1) {
						pstm2.setString(1, x);
						pstm2.setString(2, y);
						String[] tmpz = z.split("#");					
						pstm2.setString(3, tmpz[tmpz.length -1]);
						pstm2.addBatch();
						i++;
					}
					if (i % 100000 == 0) {
						System.out.println(i);						
						pstm0.executeBatch();
						pstm1.executeBatch();
						pstm2.executeBatch();
						conn.commit();
						pstm0.clearBatch();
						pstm1.clearBatch();
						pstm2.clearBatch();
					}
				}							
			}

			pstm0.executeBatch();
			pstm1.executeBatch();
			pstm2.executeBatch();
			conn.commit();
			pstm0.clearBatch();
			pstm1.clearBatch();
			pstm2.clearBatch();
			System.out.println("Database created successfully...");
			pstm0.close();
			pstm1.close();
			pstm2.close();
			stmt.close();
			conn.close();
		} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		} finally {
			// finally block used to close resources
			try {
				if (stmt != null)
					stmt.close();
			} catch (SQLException se2) {
			} // nothing we can do
			try {
				if (conn != null)
					conn.close();
			} catch (SQLException se) {
				se.printStackTrace();
			}
			
			// end finally try
		} // end try
		System.out.println("Goodbye!");
	}	
	
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
	
	static String monfile = "F:\\data\\mondial\\mondial.rdf";
	Map<String, Integer> hashNodeTable;
	ArrayList<String> reHashTable;
	ArrayList<Set<Edge>> graph;
	boolean[] visited;	
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
	

	void writeToOut(String fileName, String outfile) throws IOException {		
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();		
		graph = new ArrayList<Set<Edge>>();
		
		Model model = ModelFactory.createDefaultModel();
		InputStream in = FileManager.get().open(fileName);
		model.read(in, null, "RDF/XML-ABBREV");		
		String x, y, z;
		ResIterator subjects = model.listSubjects();		
		while (subjects.hasNext()) {
			Resource subject = subjects.next();
			// get all triples
			StmtIterator properties = subject.listProperties();
			while (properties.hasNext()) {
				Statement stmtS = properties.nextStatement();					
				RDFNode object = stmtS.getObject();
				x = subject.toString();
				y = stmtS.getPredicate().toString(); // get the predicate
				z = object.toString();
				if (object instanceof Resource) 
				{						
					int ux = tryAdd(x);
					int uz = tryAdd(z);
					graph.get(ux).add(new Edge(y, uz));
					graph.get(uz).add(new Edge(y, ux));						
				}		
			}							
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
						fopinfo.println(name[i] + " " + name[it.v]);
					}
			}
		fopinfo.close();
	}
	
	void dbdealmo() throws IOException {
		
		//writeToOut(monfile, "E://data//mondial.gg");
		/*
		String database = "mondial";
		String[] table = new String[3];
		table[0] = "graph";		
		table[1] = "keyword";
		table[2] = "type";	
		writeAll(monfile, database, table);
	*/
		WeightAssign.weightDBGen("mondial");
	}
	
	public static void main(String[] args) throws IOException {
		dbwritemondial link = new dbwritemondial();
		link.dbdealmo();
		//writeToOut();
	}
}
