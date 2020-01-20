package BANKSII;

public class MaxHeap {
	
	HeapNode []heap;	//start from 1
	int []heapLoc;
	public int heapSize;	//heap[1..heapSize] contains node in Q
		
	int keyNum;
		
	void heapInit(int node, int key) {
		int nodeNum = node + 1;
		heap = new HeapNode[nodeNum];
		heapLoc = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			heapLoc[i] = -1;		
		heapSize = 0;
		
		keyNum = key;
		
	}
	
	//swap two heap element
	void swap(int i, int j) {
		HeapNode tmp = heap[i];
		heap[i] = heap[j];
		heap[j] = tmp;
	}
	
	void push(int v, double activation) {
		heapSize++;
		heap[heapSize] = new HeapNode(v, activation);
		heapLoc[v] = heapSize;
		upModify(heapSize);
	}
		
	HeapNode pop() {
		HeapNode ans = heap[1];
		heapLoc[heap[1].v] = 0;	//0 means heap[1].v is in X
		heapLoc[heap[heapSize].v] = 1;		 
		swap(1, heapSize);
		heapSize--;
		downModify(1);
		
		return ans;
	}
	
	double topatv() {
		return heap[1].activation;
	}
	
	boolean isempty() {
		if (heapSize == 0)
			return true;
		return false;
	}
	
	boolean inQ(int v) {		
		if (heapLoc[v] >= 1)
			return true;
		return false;
	}
	
	boolean inX(int v){		
		if (heapLoc[v] == 0)
			return true;
		return false;
	}	
	
	boolean notreached(int v) {
		if (heapLoc[v] == -1)
			return true;
		return false;
	}	
	
	void change(int v, double newactivation) {
		heap[heapLoc[v]].activation = newactivation;
		upModify(heapLoc[v]);
	}
	
	//try use activation to lower v
	void tryHigher(int v, double activation) {
		if (inX(v)) return;	//poped already
		if (notreached(v)) { 	//not pushed
			push(v, activation);
			return;
		}		
		if (heap[heapLoc[v]].activation < activation)
			change(v, activation);				
	}
	
	void downModify(int loc) {
		int target = loc;
		int left = loc * 2;
		int right = loc * 2 + 1;
		if (left <= heapSize && heap[left].activation > heap[target].activation)			
			target = left;
		if (right <= heapSize && heap[right].activation > heap[target].activation)			
			target = right;
		if (target != loc) {
			heapLoc[heap[target].v] = loc;
			heapLoc[heap[loc].v] = target;
			swap(target, loc);
			downModify(target);
		}
	}
	
	void upModify(int loc) {
		if (loc > 1) {
			int target = loc / 2;
			if (heap[target].activation < heap[loc].activation) {
				heapLoc[heap[target].v] = loc;
				heapLoc[heap[loc].v] = target;
				swap(target, loc);
				upModify(target);
			}
		}
	}	
}
