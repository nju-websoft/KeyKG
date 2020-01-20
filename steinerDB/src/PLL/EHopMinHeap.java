package PLL;

public class EHopMinHeap {
	EHop []heap;	//start from 1
	int []heapLoc;
	double []dist;
	int heapSize;
	
	void heapInit(int node) {
		int nodeNum = node + 1;
		heap = new EHop[nodeNum];
		heapLoc = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			heapLoc[i] = -1;
		
		//record distance
		dist = new double[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			dist[i] = -1;
		heapSize = 0;
	}
	
	//swap two heap element
	void swap(int i, int j) {
		EHop tmp = heap[i];
		heap[i] = heap[j];
		heap[j] = tmp;
	}
	
	void push(int v, double dis, int par) {
		heapSize++;
		dist[v] = dis;
		heap[heapSize] = new EHop(v, dis, par);
		heapLoc[v] = heapSize;
		upModify(heapSize);
	}
		
	EHop pop() {
		EHop ans = heap[1];
		//dist[heap[1].v] = heap[1].dis;
		heapLoc[heap[1].v] = heapSize;
		heapLoc[heap[heapSize].v] = 1;		 
		swap(1, heapSize);
		heapSize--;
		downModify(1);
		
		return ans;
	}
	
	boolean empty() {
		if (heapSize == 0)
			return true;
		return false;
	}
	
	void change(int v, double newDis, int par) {
		heap[heapLoc[v]].dis = newDis;
		dist[v] = newDis;
		heap[heapLoc[v]].par = par;
		upModify(heapLoc[v]);
	}
	
	//try use dis to lower v
	void tryLower(int v, double dis, int par) {
		if (heapLoc[v] == -1) 	//not pushed
			push(v, dis, par);
		else {
			if (dist[v] > dis)
				change(v, dis, par);
		}			
	}
	
	void downModify(int loc) {
		int target = loc;
		int left = loc * 2;
		int right = loc * 2 + 1;
		if (left <= heapSize && heap[left].dis < heap[target].dis)			
			target = left;
		if (right <= heapSize && heap[right].dis < heap[target].dis)			
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
			if (heap[target].dis > heap[loc].dis) {
				heapLoc[heap[target].v] = loc;
				heapLoc[heap[loc].v] = target;
				swap(target, loc);
				upModify(target);
			}
		}
	}
	
}
