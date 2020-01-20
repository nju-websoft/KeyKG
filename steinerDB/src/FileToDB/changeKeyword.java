package FileToDB;

import java.io.*;
import java.util.*;

import database.dbReader;
import database.dbWriter;

public class changeKeyword {

	ArrayList<ArrayList<Integer>> invTable;	//keyword inverted table
	public HashMap<String, Integer> HashKeywordTable;	//keyword hash
	public ArrayList<String> reHashKeywordTable;	//rehash Integer to keyword
	
	//read keywrod hash table
	void readKeyMap(String DBName) {
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "keyMap");				
		int []a = new int [1];
		String []b = new String[1];
		while(dh.readKeyMap(a,b)){
			b[0] = b[0].toLowerCase();			
			if (!HashKeywordTable.containsKey(b[0])) {
				HashKeywordTable.put(b[0], reHashKeywordTable.size());
				reHashKeywordTable.add(b[0]);								
			}
		}
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
			a[0] = a[0].toLowerCase();
			invTable.get(HashKeywordTable.get(a[0])).add(b[0]);
		}
	}
	
	void dbdeal(String DBName) throws IOException{
		readKeyMap(DBName);
		readInvTable(DBName);
		//write keyword map
		dbWriter dkm = new dbWriter();
		dkm.dbInit(DBName, "keyMap");
		for (int i = 0; i < reHashKeywordTable.size(); i++)		
			dkm.insertKeyMap(i, reHashKeywordTable.get(i).toLowerCase());		
		dkm.dbWriterEnd();
		System.out.println("Keyword Map Generate end!");
		/*
		String que = "F:\\data\\linkedmdb\\query_src.txt";
		PrintWriter fop = new PrintWriter(new File(que));
		for (String st : reHashKeywordTable)
			fop.println(st.toLowerCase());
		fop.close();
		*/
		dbWriter di = new dbWriter();
		di.dbInit(DBName, "invertedTable");
		for (int i = 0; i < reHashKeywordTable.size(); i++)
			for (int j : invTable.get(i))
			di.insertInvertedTable(reHashKeywordTable.get(i).toLowerCase(), j);
		di.dbWriterEnd();
		
		System.out.println("Inverted Table Generate end!");	
	}
	
	public static void main(String[] args) throws IOException {
		changeKeyword link = new changeKeyword();
		//link.dbdeal("dbpedia201610");
		link.dbdeal("linkedmdb");
		//link.dbdeal("mondial");
	}
}
