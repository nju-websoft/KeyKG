package TestCode;

import java.io.*;
import java.util.*;

import BANKSII.GraphInfo;
import DOGST.CommonStruct;
import DOGST.DiskCommonStruct;
import DOGST.OneStarDO;
import DOGST.DOGST;
import DOGST.DOGST_D;
import BANKSII.BiSearch;
import prunedDP.GraphDP;
import prunedDP.PrunedBiHeap;
import DPBF.GraphDPBF;
import javafx.util.Pair;
import DPBF.DPBF_1;

public class TestAllGiven {
	
	String[] filest;
	ArrayList<ArrayList<String>> words = new ArrayList<ArrayList<String>>();
	void inputquery(String fileName) throws IOException {
		Scanner input1 = new Scanner(new File(fileName));		
		while (input1.hasNext()) {
			String line = input1.nextLine();
			String[] newwords = line.split(" ");
			ArrayList<String> newquery = new ArrayList<String>(); 
			for (int i = 0; i < newwords.length; i++) newquery.add(newwords[i]);
			words.add(newquery);
		}
		input1.close();
	}
	
	void outputans(String st, String fileName, ArrayList<Double> weight, ArrayList<Double> time) throws IOException {
		PrintWriter fop;
		File file = new File("E:\\data\\"+ st + "." + fileName + ".txt");
		if (!file.exists()) file.createNewFile();
		fop = new PrintWriter(file);
		
		double ttime = 0, tweight = 0;
		for (int i= 0; i < weight.size(); i++) {
			fop.println("time: "+ time.get(i) + " ms");
			fop.println("weight: "+ weight.get(i));
			ttime += time.get(i);
			tweight += weight.get(i);
		}		
		fop.println("avg time: "+ ttime/weight.size()  + " ms");
		fop.println("avg weight: "+ tweight/weight.size());
		fop.close();
	}
	
	//median 
	void outputans(String st, String fileName, List<List<Double>> weights, List<List<Double>> times) throws IOException {
		List<Double> weight = new ArrayList<>();
		List<Double> time = new ArrayList<>();
		if (times.size() == 0) return;
		int caseNum = times.get(0).size();
		int testTime = times.size();	
		for (int i= 0; i < caseNum; i++) {
			List<Pair<Double, Double>> pp = new ArrayList<>();
			for (int j = 0; j < testTime; j++)
				pp.add(new Pair<>(times.get(j).get(i), weights.get(j).get(i)));
			pp.sort(Comparator.comparing(Pair::getKey));
			double midtime = pp.get(testTime/2).getKey();
			if (testTime > 1 && testTime % 2 == 0) {
				midtime = (pp.get(testTime/2).getKey() + pp.get(testTime/2 - 1).getKey())/2;
			}
			time.add(midtime);
			weight.add(pp.get(testTime/2).getValue());
		}
		
		PrintWriter fop;
		File file = new File("E:\\data\\"+ st + "." + fileName + ".txt");
		if (!file.exists()) file.createNewFile();
		fop = new PrintWriter(file);
		
		double ttime = 0, tweight = 0;
		for (int i= 0; i < weight.size(); i++) {
			fop.println("time: "+ time.get(i) + " ms");
			fop.println("weight: "+ weight.get(i));
			ttime += time.get(i);
			tweight += weight.get(i);
		}
		fop.println("avg time: "+ ttime/weight.size()  + " ms");
		fop.println("avg weight: "+ tweight/weight.size());
		fop.close();
	}
	

	boolean testBanks = false;
	boolean testPLL = true;
	boolean testchins = true;
	boolean testPruned = false;	
	boolean testDPBF = false;
	
	
	void testDeal() throws IOException {
		filest = new String[3];
		filest[0] = "mondial";
		filest[1] = "linkedmdb";
		filest[2] = "dbpedia201610";

		String []que = new String[4];
		que[0] = "F:\\data\\mondial\\query.txt";
		que[1] = "F:\\data\\linkedmdb\\query.txt";
		que[2] = "F:\\data\\dbpedia201610\\query.txt";
		double itimes[] = {1, 1, 1};
		int dataset = 2;
		String st = filest[dataset];
		String queryfile = que[dataset];
		
		long startTime, endTime;
		words = new ArrayList<ArrayList<String>>();
		inputquery(queryfile);
		
		if (testBanks) {
			double iTime = 1;
			int testTime = 1;
			List<List<Double>> timeList = new ArrayList<>();
			List<List<Double>> weightList = new ArrayList<>();
			for (int i = 0; i < testTime; i++) {
				timeList.add(new ArrayList<>());
				weightList.add(new ArrayList<>());
			}
			GraphInfo g1 = new GraphInfo();
			g1.Init(st);
			BiSearch b1 = new BiSearch();
			
			for (int tt = 0; tt < testTime; tt++) {
				for (int i = 0; i < words.size(); i++) {
					startTime = System.currentTimeMillis();
					for (int j = 0; j < iTime; j++)
						b1.kdealRandom(g1, words.get(i));
					endTime = System.currentTimeMillis();
					weightList.get(tt).add(b1.totalweight);
					timeList.get(tt).add((double) (endTime - startTime) / iTime);					
					System.out.println(tt + " " + i + " " + b1.totalweight + " " + (double) (endTime - startTime) / iTime);
				}
			}	
			outputans(st, "BANKS-II", weightList, timeList);
		}

		if (testPLL) {
			int testTime = 1;
			CommonStruct c3 = new CommonStruct();
			c3.Init(st, true);	//test pll
			List<List<Double>> timeList = new ArrayList<>();
			List<List<Double>> weightList = new ArrayList<>();
			for (int i = 0; i < testTime; i++) {
				timeList.add(new ArrayList<>());
				weightList.add(new ArrayList<>());
			}
			OneStarDO k1 = new OneStarDO(c3);	// KeyKG-pll
			double kgtime = itimes[dataset];
			for (int tt = 0; tt < testTime; tt++) {
				for (int i = 0; i < words.size(); i++) {
					startTime = System.currentTimeMillis();
					for (int j = 0; j < kgtime; j++)
						k1.kdealRandom(c3, words.get(i), 2);
					endTime = System.currentTimeMillis();
					weightList.get(tt).add(k1.totalweight);
					timeList.get(tt).add((double) (endTime - startTime) / kgtime);
					//System.out.println(tt + " " + i + " " + k1.totalweight + " " + (double) (endTime - startTime) / kgtime);
				}
			}
			//outputans(st, "KeyKG-PLL", weightList, timeList);
		}
		
		if (testchins) {
			
			CommonStruct c1 = new CommonStruct();
			c1.Init(st);
			
			double iTime = 1;
			int testTime = 1;
					
			boolean test_KeyKG = true;
			boolean test_KeyKGP = false;
			boolean test_D = false;
			
			if (test_KeyKG) {
				List<List<Double>> timeList = new ArrayList<>();
				List<List<Double>> weightList = new ArrayList<>();
				for (int i = 0; i < testTime; i++) {
					timeList.add(new ArrayList<>());
					weightList.add(new ArrayList<>());
				}
				OneStarDO k1 = new OneStarDO(c1);	// KeyKG
				double kgtime = itimes[dataset];
				for (int tt = 0; tt < testTime; tt++) {
					for (int i = 0; i < words.size(); i++) {
						startTime = System.currentTimeMillis();
						for (int j = 0; j < kgtime; j++)
							k1.kdealRandom(c1, words.get(i), 2);
						endTime = System.currentTimeMillis();
						weightList.get(tt).add(k1.totalweight);
						timeList.get(tt).add((double) (endTime - startTime) / kgtime);
						//System.out.println(tt + " " + i + " " + k1.totalweight + " " + (double) (endTime - startTime) / kgtime);
					}
				}
				//outputans(st, "KeyKG", weightList, timeList);
			}
			
			if (test_KeyKGP) {
				List<List<Double>> timeList = new ArrayList<>();
				List<List<Double>> weightList = new ArrayList<>();
				for (int i = 0; i < testTime; i++) {
					timeList.add(new ArrayList<>());
					weightList.add(new ArrayList<>());
				}
				DOGST k3 = new DOGST(c1); // keyKG+
				
				for (int tt = 0; tt < testTime; tt++) {
					for (int i = 0; i < words.size(); i++) {
						startTime = System.currentTimeMillis();
						for (int j = 0; j < iTime; j++)
							k3.kdealRandom(c1, words.get(i), 2);
						endTime = System.currentTimeMillis();
						weightList.get(tt).add(k3.totalweight);
						timeList.get(tt).add((double) (endTime - startTime) / iTime);
						System.out.println(tt + " " + i + " " + k3.totalweight + " " + (double) (endTime - startTime) / iTime);
					}
				}
				outputans(st, "KeyKG+", weightList, timeList);
			}			
		
			if (test_D) {
				List<List<Double>> timeList = new ArrayList<>();
				List<List<Double>> weightList = new ArrayList<>();
				for (int i = 0; i < testTime; i++) {
					timeList.add(new ArrayList<>());
					weightList.add(new ArrayList<>());
				}				
				DiskCommonStruct c2 = new DiskCommonStruct();
				c2.Init(st);
				DOGST_D k1 = new DOGST_D(c2); // KeyKGD
				double disktime = itimes[dataset];
				for (int tt = 0; tt < testTime; tt++) {
					for (int i = 0; i < words.size(); i++) {
						startTime = System.currentTimeMillis();
						for (int j = 0; j < disktime; j++)
							k1.kdealRandom(c2, words.get(i), 2);
						endTime = System.currentTimeMillis();
						weightList.get(tt).add(k1.totalweight);
						timeList.get(tt).add((double) (endTime - startTime) / disktime);
						System.out.println(tt + " " + i + " " + k1.totalweight + " " + (double) (endTime - startTime) / disktime);
					}
				}
				outputans(st, "KeyKG+-D", weightList, timeList);
			}
			
			
		}
		
		if (testPruned) {
			double iTime = 1;
			int testTime = 1;
			List<List<Double>> timeList = new ArrayList<>();
			List<List<Double>> weightList = new ArrayList<>();
			for (int i = 0; i < testTime; i++) {
				timeList.add(new ArrayList<>());
				weightList.add(new ArrayList<>());
			}
						
			GraphDP g1 = new GraphDP();
			g1.Init(st);
			PrunedBiHeap p1 = new PrunedBiHeap();
			Double ap = 1.0;
			for (int tt = 0; tt < testTime; tt++) {
				for (int i = 0; i < words.size(); i++) {		
					g1.givenQueryWord(words.get(i));				
					startTime = System.currentTimeMillis();				
					for (int j = 0; j < iTime; j++)
						p1.kdeal(g1, ap, 3600*1000);				
					endTime = System.currentTimeMillis();
					weightList.get(tt).add(p1.best);
					timeList.get(tt).add((double)(endTime - startTime)/iTime);
					System.out.println(tt + " " + i + " " + p1.best + " " + ((double)(endTime - startTime)/1000));					
				}
			}
			outputans(st, "PrunedDP++", weightList, timeList);
		}
		
		if (testDPBF) {
			double iTime = 1;
			int testTime = 1;
			List<List<Double>> timeList = new ArrayList<>();
			List<List<Double>> weightList = new ArrayList<>();
			for (int i = 0; i < testTime; i++) {
				timeList.add(new ArrayList<>());
				weightList.add(new ArrayList<>());
			}
			GraphDPBF g1 = new GraphDPBF();
			g1.Init(st);
			DPBF_1 d1 = new DPBF_1();
			for (int tt = 0; tt < testTime; tt++) {
				for (int i = 0; i < words.size(); i++) {
					g1.givenQueryWord(words.get(i));
					startTime = System.currentTimeMillis();
					for (int j = 0; j < iTime; j++)
						d1.kdeal(g1);
					endTime = System.currentTimeMillis();		
					weightList.get(tt).add(d1.best);
					timeList.get(tt).add((double)(endTime - startTime)/iTime);
				}			
			}
			outputans(st, "DPBF", weightList, timeList);
		}
		
		
	}

	public static void main(String[] args) throws IOException {
		TestAllGiven t1 = new TestAllGiven();
		t1.testDeal();
		// t1.testDealDisk();
	}
}
