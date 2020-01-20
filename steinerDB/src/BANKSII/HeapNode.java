package BANKSII;

public class HeapNode {
	int v;
	double activation;		
	HeapNode(int v, double activation){
		this.v = v;
		this.activation = activation;		
	}
	HeapNode(int v, double activation, int depth){
		this.v = v;
		this.activation = activation;		
	}
}
