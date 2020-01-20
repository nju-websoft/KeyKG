#ifndef __CentInt__
#define __CentInt__
#include "pllInt.h"
//#include "pllThread.h"
#include "UsedSt.h"
#include <random>
#include <algorithm>
using namespace std;

typedef struct valcentInt {
	int name;
	double f;
	valcentInt(int name, double f) {
		this->name = name;
		this->f = f;
	}
} valcentInt;

static bool lessSortcentInt(valcentInt v1, valcentInt v2) {
	if (fabs(v1.f - v2.f) > 1e-8)
		return v1.f > v2.f;
	else
		return v1.name > v2.name;
}
class centInt {
public:
	centInt();
	void deal(string, bool, int);
	void chlV1(int);
	void chlV2(int);

private:
	string file;
	int nodeNum;
	pllInt p1;

	vector<valcentInt> valList;
	vector<int> oldName;
	vector<int> newName;
};
#endif
