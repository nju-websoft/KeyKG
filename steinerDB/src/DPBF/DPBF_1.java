package DPBF;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.TreeMap;


public class DPBF_1 {
	GraphDPBF g1;
	int realName[];
	int nodeNum;
	int keyNum;	
	BitSet P;
	MST m1;
	PriorityQueue<CandiTree> QT;
	Map<DTree, DTreeInfo> DT;
	public double best = 0;
	
	int getloc(BitSet X) {
		int ans = 0;
		for (int i = X.nextSetBit(0); i >= 0; i = X.nextSetBit(i + 1))
		      ans += X.get(i) ? (1 << i) : 0;
		return ans;
	}
	
	
	double DPBFDeal() throws FileNotFoundException {		
		//record P in the form of BitSet 
		P = new BitSet(keyNum);
		P.flip(0, keyNum);
		
		realName = new int[keyNum];
		for (int i = 0; i < keyNum; i++)
			realName[i] = g1.HashKeywordTable.get(g1.keywordList.get(i));
		double keynode = 0;
		for (int i = 0; i < keyNum; i++)	//record all 
			keynode = keynode + g1.invTable.get(realName[i]).size();		
		
		m1 = new MST(nodeNum);
		QT = new PriorityQueue<CandiTree>();	
		DT = new TreeMap<DTree, DTreeInfo>();
		
		PrintWriter output = new PrintWriter(new File("E:\\data\\test.txt"));
		//line 3-5
		for (int i = 0; i < keyNum; i++)
			for (int it : g1.invTable.get(realName[i]))
				QT.add(new CandiTree(it, keyNum, i));
		while (!QT.isEmpty()) {	//line 6
			CandiTree ct = QT.poll();			
			if (DT.containsKey(new DTree(ct.v, ct.X)))
				continue;
			DT.put(new DTree(ct.v, ct.X), new DTreeInfo(ct.weight, ct.edges));
			output.println(ct.v + " " + getloc(ct.X)  + " " + ct.weight);
			//line 7
			if (ct.X.equals(P)) {
				best = ct.weight;
				System.out.println(best);
				//for(EdgeDPBF it : ct.edges) System.out.println(it.u + " " + it.v + " " + it.weight);
				output.close();
				return keynode;				
			}
			int v = ct.v;
			for (EdgeDPBF it : g1.graph.get(v)) {  	//line 9
				DTreeInfo dti = m1.prim(edgesMerge(ct.edges, it));
				QT.add(new CandiTree(it.v, ct.X, dti.cost, dti.edges));	//line 10
			}
				
			
			BitSet Xbar = (BitSet)P.clone();			
			Xbar.xor(ct.X);
			findSubset(ct, ct.X, Xbar, 0);
		}
		return keynode;
	}
	
	//line 13-17
	void findSubset(CandiTree ct, BitSet X, BitSet Xbar, int loc) {
		if (Xbar.cardinality() > 0) {
			DTree newdt = new DTree(ct.v, Xbar);
			if (DT.containsKey(newdt)) {
				X.xor(Xbar);
				DTree candt = new DTree(ct.v, X);
				if (!DT.containsKey(candt)) {
					DTreeInfo dti = m1.prim(edgesMerge(ct.edges, DT.get(newdt).edges));
					QT.add(new CandiTree(ct.v, X, dti.cost, dti.edges));					
				}
				X.xor(Xbar);
			}
		}
		for (int i = Xbar.nextSetBit(loc); i >= 0; i = Xbar.nextSetBit(i + 1)) {
			Xbar.flip(i);
			findSubset(ct, X, Xbar, i + 1);
			Xbar.flip(i);
		}
	}	
	
	ArrayList<EdgeDPBF> edgesMerge(ArrayList<EdgeDPBF> edges1, EdgeDPBF edge2) {
		ArrayList<EdgeDPBF> edges = new ArrayList<EdgeDPBF>();
		for (EdgeDPBF it : edges1) edges.add(it);
		edges.add(edge2);
		return edges;
	}
	
	ArrayList<EdgeDPBF> edgesMerge(ArrayList<EdgeDPBF> edges1, ArrayList<EdgeDPBF> edges2) {
		ArrayList<EdgeDPBF> edges = new ArrayList<EdgeDPBF>();
		for (EdgeDPBF it : edges1) edges.add(it);
		for (EdgeDPBF it : edges2) edges.add(it);		
		return edges;
	}
	
	public double kdeal(GraphDPBF gg) throws FileNotFoundException {
		g1 = gg;
		nodeNum = g1.nodeNum;
		keyNum = g1.keyNum;		
		return DPBFDeal();
	}
}
