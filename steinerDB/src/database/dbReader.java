package database;

import java.sql.*;
import java.util.ArrayList;

import PLL.EHop;

public class dbReader {
	static final String JDBC_DRIVER = "com.mysql.jdbc.Driver";
	static final String DB_URL = "jdbc:mysql://localhost:3306/?user=root&characterEncoding=utf8&useCursorFetch=true";
	static final String USER = "root";
	static final String PASS = "shiyuxuan";
	String DBName;
	String tableName;

	Connection conn = null;
	java.sql.Statement stmt = null;
	ResultSet rs = null;
	ModeEnum mode = ModeEnum.Undefined;
	public void dbInit(String database, String table) {
		DBName = database;
		tableName = table;

		String sql;

		try {
			Class.forName("com.mysql.jdbc.Driver");
			conn = DriverManager.getConnection(DB_URL, USER, PASS);

			stmt = conn.createStatement();			
			
			sql = "CREATE DATABASE IF NOT EXISTS " + database;			
			stmt.executeUpdate(sql);
			
			// use database
			sql = "USE " + database;
			stmt.executeUpdate(sql);

			sql = "SELECT * from " + table;
			stmt.setFetchSize(100000);
			rs = stmt.executeQuery(sql);
			if (rs.wasNull())
				System.out.println("Open table " + table + " fail");
			
			if (table.equals("subName")||table.equals("subNamePLL"))
				mode = ModeEnum.SubName;
			if (table.equals("graph"))
				mode = ModeEnum.Graph;
			if (table.equals("hubLabel")||table.equals("hubLabelPLL"))
				mode = ModeEnum.HubLabel;
			if (table.equals("invertedTable")||table.equals("invertedTablePLL"))
				mode = ModeEnum.InvertedTable;
			if (table.equals("edgeWeight"))
				mode = ModeEnum.EdgeWeight;
			if (table.equals("keyword"))
				mode = ModeEnum.Keyword;
			if (table.equals("keyMap"))
				mode = ModeEnum.KeyMap;			
			
		} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}

	}

	/*
	public void dbInit(String database, String table, int u) {
		DBName = database;
		tableName = table;
		String sql;
		try {
			Class.forName("com.mysql.jdbc.Driver");
			conn = DriverManager.getConnection(DB_URL, USER, PASS);
			stmt = conn.createStatement();						
			sql = "CREATE DATABASE IF NOT EXISTS " + database;			
			stmt.executeUpdate(sql);			
			// use database
			sql = "USE " + database;
			stmt.executeUpdate(sql);
			sql = "SELECT * from " + table + " WHERE u = "+ u;			
			stmt.setFetchSize(100000);
			rs = stmt.executeQuery(sql);
			if (rs.wasNull())
				System.out.println("fetch table " + table + " fail");

			if (table.equals("hubLabel"))
				mode = ModeEnum.HubLabel;			
		} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}
	}
	*/
	
	public void dbInitDisk(String database, String table) {
		DBName = database;
		tableName = table;
		String sql;
		try {
			Class.forName("com.mysql.jdbc.Driver");
			conn = DriverManager.getConnection(DB_URL, USER, PASS);
			stmt = conn.createStatement();						
			sql = "CREATE DATABASE IF NOT EXISTS " + database;			
			stmt.executeUpdate(sql);			
			// use database
			sql = "USE " + database;
			stmt.executeUpdate(sql);
			
			} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}
	}
	
	public void dbGet(int u, ArrayList<EHop> newla) {
		String sql;
		try {					
			sql = "SELECT * from  hubLabel WHERE u = "+ u;
			stmt.setFetchSize(500);
			rs = stmt.executeQuery(sql);
			if (rs.wasNull())
				System.out.println("fetch table hubLabel fail");
			while (rs.next())
				newla.add(new EHop(rs.getInt("v"), rs.getDouble("dis"), rs.getInt("par")));
			rs.close();
			} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}
		
	}
	
	public void dbGet(int u, int v,int[] singleP, double[] single) {
		String sql;
		try {
					
			sql = "SELECT * from  hubLabel WHERE u = "+ u + " AND v = " + v;
			stmt.setFetchSize(500);
			rs = stmt.executeQuery(sql);
			if (rs.wasNull())
				System.out.println("fetch table hubLabel fail");
			while (rs.next()) {
				single[0] = rs.getDouble("dis");
				singleP[0] = rs.getInt("par");
				//System.out.println(rs.getInt("u")+" "+ rs.getInt("v")+" " + single[0]);
			}
			rs.close();
			} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}
		
	}
	
	public boolean readGraph(String triples[]) {
		if (mode != ModeEnum.Graph)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				triples[0] = rs.getString("subject");
				triples[1] = rs.getString("predicate");
				triples[2] = rs.getString("object");
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}

	public boolean readHubLabel(int[] triples, double[] single) {
		if (mode != ModeEnum.HubLabel)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				
				triples[0] = rs.getInt("u");
				triples[1] = rs.getInt("v");
				single[0] = rs.getDouble("dis");
				triples[2] = rs.getInt("par");				
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean readInvertedTable(String[] singleK, int[] singleN) {
		if (mode != ModeEnum.InvertedTable)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				
				singleK[0] = rs.getString("keyword");
				singleN[0] = rs.getInt("contained");	
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean readEdgeWeight(double[] single) {
		if (mode != ModeEnum.EdgeWeight)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line				
				single[0] = rs.getDouble("weight");				
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean readKeyword(String triples[]) {
		if (mode != ModeEnum.Keyword)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				triples[0] = rs.getString("subject");
				triples[1] = rs.getString("predicate");
				triples[2] = rs.getString("object");
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean readSubName(String[] singleK, int[] singleN) {
		if (mode != ModeEnum.SubName)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				
				singleK[0] = rs.getString("subject");
				singleN[0] = rs.getInt("newName");	
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean readKeyMap(int[] singleN, String[] singleK) {		
		if (mode != ModeEnum.KeyMap)
			return false;
		try {
			if (rs.isClosed())
				return false;
			if (rs.next()) { // get a new line
				
				singleN[0] = rs.getInt("node");
				singleK[0] = rs.getString("keyword");				
				return true;
			} else {
				rs.close();
				stmt.close();
				conn.close();
				return false;
			}
		} catch (SQLException e) {			
			e.printStackTrace();
		}
		return false;
	}
	
	public static void main(String[] args) {

	}
}
