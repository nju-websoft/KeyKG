package DPBF;

import java.util.ArrayList;

public class DTreeInfo {
	double cost;
	ArrayList<EdgeDPBF> edges;
	DTreeInfo(double cost, ArrayList<EdgeDPBF> edges){
		this.cost = cost;
		this.edges = new ArrayList<EdgeDPBF>();
		if (edges == null)
			return;
		for (EdgeDPBF it : edges)
			this.edges.add(new EdgeDPBF(it));
	}
}
