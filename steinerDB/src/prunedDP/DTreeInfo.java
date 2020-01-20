package prunedDP;

import java.util.ArrayList;

public class DTreeInfo {
	double cost;
	ArrayList<EdgeDP> edges;
	DTreeInfo(double cost, ArrayList<EdgeDP> edges){
		this.cost = cost;
		this.edges = new ArrayList<EdgeDP>();
		if (edges == null)
			return;
		for (EdgeDP it : edges)
			this.edges.add(new EdgeDP(it));
	}
	
	void print(){
		System.out.println(cost);
		for (EdgeDP it : edges) {
			System.out.println(it.u + " " + it.v + " " + it.weight);
		}
	}
}
