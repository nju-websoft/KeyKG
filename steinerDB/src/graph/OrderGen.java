package graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class OrderGen {
	public void orderGen1(WeightedGraph w1, int[] newName, int[] oldName) {
		ArrayList<Hop> nodeDegree = new ArrayList<Hop>(); 
		for (int i = 0; i < w1.nodeNum; i++)
			nodeDegree.add(new Hop(i, w1.graph.get(i).size()));
		Collections.sort(nodeDegree, new HopCompare());		
		for (int i = 0; i < w1.nodeNum; i++)
			oldName[i] = nodeDegree.get(i).v;
		for (int i = 0; i < w1.nodeNum; i++)
			newName[oldName[i]] = i;
	}

	public void orderGen2(WeightedGraph w1, int[] newName, int[] oldName) {
		ArrayList<Hop> nodeDegree = new ArrayList<Hop>(); 
		for (int i = 0; i < w1.nodeNum; i++)
			nodeDegree.add(new Hop(i, w1.graph.get(i).size()));
		Collections.sort(nodeDegree, new HopCompare());
		//use top 100 node to calculate betweeness
		Brandes b1 = new Brandes(w1);
		int itTime = 200;
		if (itTime > w1.nodeNum)
			itTime = w1.nodeNum;
		
		for (int i = 0; i < itTime; i++) {
			b1.dfs(nodeDegree.get(i).v);
			if (i % 20 == 0)
				System.out.println("bc " + i +" end");
		}
		
		//sort node according to node bc
		ArrayList<Hop> nodeCen = new ArrayList<Hop>(); 
		for (int i = 0; i < w1.nodeNum; i++)
			nodeCen.add(new Hop(i,b1.centrality[i]));
		Collections.sort(nodeCen, new HopCompare());	
		
		for (int i = 0; i < w1.nodeNum; i++)
			oldName[i] = nodeCen.get(i).v;
		for (int i = 0; i < w1.nodeNum; i++)
			newName[oldName[i]] = i;
	}
	
	public void orderGenRand(WeightedGraph w1, int[] newName, int[] oldName) {		
		ArrayList<Integer> solution = new ArrayList<>();
		for (int i = 0; i < w1.nodeNum; i++) 
		    solution.add(i);
		Collections.shuffle(solution);
		//use top 100 node to calculate betweeness
		Brandes b1 = new Brandes(w1);
		int itTime = 200;
		if (itTime > w1.nodeNum)
			itTime = w1.nodeNum;
		for (int i = 0; i < itTime; i++) {
			b1.dfs(solution.get(i));
			if (i % 20 == 0)
				System.out.println("bc " + i +" end");
		}
		
		//sort node according to node bc
		ArrayList<Hop> nodeCen = new ArrayList<Hop>(); 
		for (int i = 0; i < w1.nodeNum; i++)
			nodeCen.add(new Hop(i,b1.centrality[i]));
		Collections.sort(nodeCen, new HopCompare());
		for (int i = 0; i < w1.nodeNum; i++)
			oldName[i] = nodeCen.get(i).v;
		for (int i = 0; i < w1.nodeNum; i++)
			newName[oldName[i]] = i;
	}
	
	public static void main(String[] args) throws IOException {
	}
}
