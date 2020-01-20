#ifndef __CENT__
#define __CENT__
#include "pll.h"
#include "UsedSt.h"
#include <random>
using namespace std;

typedef struct valCent {
	int name;
	double f;
	valCent(int name, double f) {
		this->name = name;
		this->f = f;
	}
} valCent;

static bool lessSortCent(valCent v1, valCent v2) {
	if (fabs(v1.f - v2.f) > 1e-8)
		return v1.f > v2.f;
	else
		return v1.name > v2.name;
}
class cent {
public:
	cent();
	void deal(string, bool, int);
	void deal(string, string, int);
	void chlV1(int);
	void chlV2(int);

private:
	string file;
	int nodeNum;
	pll p1;
	
	vector<valCent> valList;	
	vector<int> oldName;
	vector<int> newName;
};
#endif
