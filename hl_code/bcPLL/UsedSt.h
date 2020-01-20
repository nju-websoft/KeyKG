#ifndef __USEDST__
#define __USEDST__
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cstring>
#include <set>
#include <time.h> 
#include <math.h>

enum Mode{ RUN_PLLD = 0, RUN_PLLR, GENKG, GENRD, RUN_PLLPAL};
extern  Mode CodeMode;

typedef struct hop {
	int goal;
	double len;
	hop(int goal, double len) {
		this->goal = goal;
		this->len = len;
	}
	hop(){
		this->goal = -1;
		this->len = -1;
	}
	friend bool operator < (struct hop const &a, struct hop const &b)
	{
		if (a.goal == b.goal){
			return a.len < b.len;
		}
		else{
			return a.goal < b.goal;
		}
	}
} hop;

typedef struct hopInt {
    int goal;
    int len;
    hopInt(int goal, int len) {
        this->goal = goal;
        this->len = len;
    }
    hopInt(){
        this->goal = -1;
        this->len = -1;
    }
    friend bool operator < (struct hopInt const &a, struct hopInt const &b)
    {
        if (a.goal == b.goal)
        {
            return a.len < b.len;
        }
        else
        {
            return a.goal < b.goal;
        }
    }
} hopInt;

typedef struct hoppath {
	int goal;
	double len;
	int father;
	hoppath(int goal, double len, int father) {
		this->goal = goal;
		this->len = len;
	}
	hoppath(){
		this->goal = -1;
		this->len = -1;
	}
	bool operator < (struct hoppath const &b){
		return father < b.father;
	}
} hoppath;
#endif
