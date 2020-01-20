package database;

import java.sql.*;

public class dbWriter {
	static final String JDBC_DRIVER = "com.mysql.jdbc.Driver";
	static final String DB_URL = "jdbc:mysql://localhost:3306/?user=root&characterEncoding=utf8";
	static final String USER = "root";
	static final String PASS = "shiyuxuan";
	String DBName;
	String tableName;

	Connection conn = null;
	java.sql.Statement stmt = null;	
	PreparedStatement pstm = null;
	int count;
	ModeEnum mode;
	int breakNum = 200000;
	public void dbInit(String database, String table) {
		DBName = database;
		tableName = table;

		String sql;
		count = 0;
		mode = ModeEnum.Undefined;
		try {
			Class.forName("com.mysql.jdbc.Driver");			
			conn = DriverManager.getConnection(DB_URL, USER, PASS);

			stmt = conn.createStatement();
			
			sql = "CREATE DATABASE IF NOT EXISTS " + database;			
			stmt.executeUpdate(sql);
			// use database
			sql = "USE " + database;
			stmt.executeUpdate(sql);

			sql = "DROP TABLE IF EXISTS "+ table;
			stmt.executeUpdate(sql);
			
			if (table.equals("subName") || table.equals("subNamePLL")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(subject VARCHAR(1023),"			
					+ " newName INTEGER)";
				stmt.executeUpdate(sql);						
				sql = "INSERT INTO "+table+"(subject, newName) VALUES (?,?)";
				mode = ModeEnum.SubName;
			}
			
			if (table.equals("hubLabel") || table.equals("hubLabelPLL")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(u INTEGER,"
					+ " v INTEGER, dis DOUBLE, par INTEGER)";
				stmt.executeUpdate(sql);						
				sql = "INSERT INTO "+table+"(u, v, dis, par) VALUES (?,?,?,?)";
				mode = ModeEnum.HubLabel;
			}
			
			if (table.equals("invertedTable")||table.equals("invertedTablePLL")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(keyword VARCHAR(1023),"			
						+ " contained INTEGER)";
				stmt.executeUpdate(sql);						
				sql = "INSERT INTO "+table+"(keyword, contained) VALUES (?,?)";
				mode = ModeEnum.InvertedTable;
			}
			
			if (table.equals("edgeWeight")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(weight DOUBLE)";
				stmt.executeUpdate(sql);						
				sql = "INSERT INTO "+table+"(weight) VALUES (?)";
				mode = ModeEnum.EdgeWeight;
			}
			
			if (table.equals("graph")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(subject VARCHAR(1023),"
						+ "predicate VARCHAR(1023),"
						+ " object VARCHAR(1023))";
				stmt.executeUpdate(sql);				
				sql = "INSERT INTO "+table+"(subject, predicate, object) VALUES (?,?,?)";
				mode = ModeEnum.Graph;
			}
			
			if (table.equals("keyword")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(subject VARCHAR(1023),"
						+ "predicate VARCHAR(1023),"
						+ " object VARCHAR(1023))";
				stmt.executeUpdate(sql);				
				sql = "INSERT INTO "+table+"(subject, predicate, object) VALUES (?,?,?)";
				mode = ModeEnum.Keyword;
			}
			
			if (table.equals("keyMap")) {
				sql = "CREATE TABLE IF NOT EXISTS "+table+"(node INTEGER,"
						+ " keyword VARCHAR(1023))";
				stmt.executeUpdate(sql);						
				sql = "INSERT INTO "+table+"(node, keyword) VALUES (?,?)";
				mode = ModeEnum.KeyMap;
			}
			
			pstm = conn.prepareStatement(sql);
			conn.setAutoCommit(false);
			
				
		} catch (SQLException se) {
			// Handle errors for JDBC
			se.printStackTrace();
		} catch (Exception e) {
			// Handle errors for Class.forName
			e.printStackTrace();
		}

	}

	//inset subName into database
	public void insertSubName(String subject, int newName) {
		if (mode != ModeEnum.SubName)
			return;
		count++;
		try {
			pstm.setString(1, subject);
			pstm.setInt(2, newName);
			pstm.addBatch();
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {			
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	//inset hubLabel into database
	public void insertHubLabel(int u, int v, double dis, int par) {
		if (mode != ModeEnum.HubLabel)
			return;
		count++;
		try {
			pstm.setInt(1, u);
			pstm.setInt(2, v);
			pstm.setDouble(3, dis);
			pstm.setInt(4, par);
			pstm.addBatch();
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {			
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//inset invertTabel into database
	public void insertInvertedTable(String keyword, int contained) {
		if (mode != ModeEnum.InvertedTable)
			return;
		count++;
		try {
			pstm.setString(1, keyword);
			pstm.setInt(2, contained);
			pstm.addBatch();
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	//inset edgeWeight into database
	public void insertEdgeWeight(double weight) {
		if (mode != ModeEnum.EdgeWeight)
			return;
		count++;
		try {
			pstm.setDouble(1, weight);
			pstm.addBatch();
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//inset graph into database
	public void insertGraph(String x, String y, String z) {
		if (mode != ModeEnum.Graph)
			return;
		count++;
		try {
			pstm.setString(1, x);
			pstm.setString(2, y);
			pstm.setString(3, z);
			pstm.addBatch();			
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void insertKeyword(String x, String y, String z) {
		if (mode != ModeEnum.Keyword)
			return;
		count++;
		try {
			pstm.setString(1, x);
			pstm.setString(2, y);
			pstm.setString(3, z);
			pstm.addBatch();			
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void insertKeyMap(int node, String keyword) {
		if (mode != ModeEnum.KeyMap)
			return;
		count++;
		try {
			pstm.setInt(1, node);
			pstm.setString(2, keyword);			
			pstm.addBatch();
			if (count % breakNum == 0) {
				System.out.println(count);
				pstm.executeBatch();
				conn.commit();
				pstm.clearBatch();
			}
		} catch (SQLException se) {			
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//write all 
	public void dbWriterEnd() {
		try {
			pstm.executeBatch();
			conn.commit();
			pstm.clearBatch();
			
			if (mode == ModeEnum.HubLabel && tableName.equals("hubLabel"))		
				stmt.executeUpdate("CREATE INDEX uv ON "+ tableName +" (u,v)");
			pstm.close();
			stmt.close();
			conn.close();
		} catch (SQLException se) {
			se.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {

	}
}
