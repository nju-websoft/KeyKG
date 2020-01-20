package TestCode;

import java.io.IOException;
import java.lang.management.*;
import java.util.ArrayList;
import java.util.Scanner;

import PLL.EHop;
import database.dbReader;

public class MemoryAnalyse {
	ArrayList<ArrayList<EHop>> label; // 2 hop label
	//read hub label
	void readHubLabel(String DBName){
		label = new ArrayList<ArrayList<EHop>>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "hubLabel");
		int []a = new int[3];
		double []b = new double [1];
		while(dh.readHubLabel(a,b)){
			if (label.size() <= a[0])
				label.add(new ArrayList<EHop>());
			label.get(a[0]).add(new EHop(a[1],b[0],a[2]));
		}
	}
	
	void testDeal() {
		String[] filest = new String[4];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		filest[2] = "dbpedia201610";
		
		String st = filest[1];
		MemoryMXBean bean = ManagementFactory.getMemoryMXBean();

		MemoryUsage memoryUsage = bean.getHeapMemoryUsage();
		long l = memoryUsage.getUsed();
		readHubLabel(st);

		long r = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getUsed();
		System.out.println(l/1024/1024);
		System.out.println(r/1024/1024);
		System.out.println((r-l)/1024/1024);
		
	}
	public static void main(String[] args) throws IOException {
		MemoryAnalyse t1 = new MemoryAnalyse();
		t1.testDeal();
		// t1.testDealDisk();
	}
}
