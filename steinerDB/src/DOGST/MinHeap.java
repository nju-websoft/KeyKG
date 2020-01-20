package DOGST;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import DOGST.DOGSTSC.setCoverSt;

public class MinHeap {
	
	setCoverSt []heap;	//start from 1
	int []heapLoc;
	Map<Integer, Integer> hashName;
	public int heapSize;	//heap[1..heapSize] contains node in Q
	
	MinHeap(ArrayList<setCoverSt> scsts){
		heapSize = 0;
		heap = new setCoverSt[scsts.size() + 1];
		for (int i = 0; i < scsts.size(); i++) heap[i] = null;
		hashName = new HashMap<Integer, Integer>();
		for (int i = 0; i < scsts.size(); i++)
			hashName.put(scsts.get(i).v, hashName.size());
		heapLoc = new int[scsts.size()];
		for (int i = 0; i < scsts.size(); i++) heapLoc[i] = -1;
		for (setCoverSt scst : scsts)
			push(scst);
	}
	
	//swap two heap element
	void swap(int i, int j) {
		setCoverSt tmp = heap[i];
		heap[i] = heap[j];
		heap[j] = tmp;
	}
	
	void push(setCoverSt scst) {
		heapSize++;
		heap[heapSize] = scst;
		heapLoc[hashName.get(scst.v)] = heapSize;
		upModify(heapSize);
	}
		
	setCoverSt poll() {
		if (heapSize == 0) return null;
		setCoverSt ans = heap[1];
		heapLoc[hashName.get(heap[1].v)] = 0;
		heapLoc[hashName.get(heap[heapSize].v)] = 1;		 
		swap(1, heapSize);
		heapSize--;
		downModify(1);
		
		return ans;
	}
	
	
	boolean isEmpty() {
		if (heapSize == 0)
			return true;
		return false;
	}
	
	void change(int v,int i) {
		if (heapLoc[hashName.get(v)] <= 0) return;
		setCoverSt scst = heap[heapLoc[hashName.get(v)]];
		if (!scst.X.get(i)) System.out.println("Heap change wrong!");
		scst.X.clear(i);
		if (scst.X.cardinality() == 0)
			scst.prio = Double.MAX_VALUE;
		else
			scst.prio = scst.weight / (scst.X.cardinality());
		downModify(heapLoc[hashName.get(v)]);
	}
		
	void downModify(int loc) {
		int target = loc;
		int left = loc * 2;
		int right = loc * 2 + 1;
		if (left <= heapSize && heap[left].prio < heap[target].prio)			
			target = left;
		if (right <= heapSize && heap[right].prio < heap[target].prio)			
			target = right;
		if (target != loc) {
			heapLoc[hashName.get(heap[target].v)] = loc;
			heapLoc[hashName.get(heap[loc].v)] = target;
			swap(target, loc);
			downModify(target);
		}
	}
	
	void upModify(int loc) {
		if (loc > 1) {
			int target = loc / 2;
			if (heap[target].prio > heap[loc].prio) {
				heapLoc[hashName.get(heap[target].v)] = loc;
				heapLoc[hashName.get(heap[loc].v)] = target;
				swap(target, loc);
				upModify(target);
			}
		}
	}	
}
