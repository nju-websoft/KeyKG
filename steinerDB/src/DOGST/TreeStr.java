package DOGST;

public class TreeStr {
	public double weight;
	public int root;
	public int []child;
	
	TreeStr(){
		weight = 0;
		root = -1;	
	}
	
	TreeStr(int len){
		weight = 0;
		root = -1;
		child = new int[len];
		child[0] = -10;
		for (int i = 1; i < len; i++)
			child[i] = 0;
	}
	
	TreeStr(TreeStr t1){
		weight = t1.weight;
		root = t1.root;
		child = new int[t1.child.length];
		for (int i = 0; i < t1.child.length; i++)
			child[i] = t1.child[i];
	}
}
