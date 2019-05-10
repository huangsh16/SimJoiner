#ifndef __EXP2_SIMJOINER_H__
#define __EXP2_SIMJOINER_H__

#include <vector>
#include <utility>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <ctime>
#include <fstream>
#include <unordered_map>
#include <set>
using namespace std;

const int SUCCESS = 0;
const int FAILURE = 1;
const int STRSIZE = 260;
const int LEN = 10;
const int INF = (int)1e9;
const int STRNUM = 1000000;
#define pb push_back
#define A first
#define B second
#define mp make_pair
typedef unsigned long long ull;

template <typename IDType, typename SimType>
struct JoinResult {
    IDType id1;
    IDType id2;
    SimType s;
    bool operator < (const JoinResult &jr) const {
	    return id1 < jr.id1 || (id1 == jr.id1 && id2 < jr.id2);
	}
	bool operator == (const JoinResult &jr) const {
		return id1 == jr.id1 && id2 == jr.id2;
	}
	JoinResult(){}
	JoinResult(IDType a, IDType b, SimType c) : id1(a), id2(b), s(c) {}
};

typedef JoinResult<unsigned, double> JaccardJoinResult;
typedef JoinResult<unsigned, unsigned> EDJoinResult;

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();

    int joinTimes;
	int *isAppear;
	int idNum;
//jacc
	unordered_map<ull, int> strHashToId;
	int *idToCnt;
    vector<vector<int> > idVectorVector[2];
    vector<vector<pair<int, int>> > TmpVec;
    vector<int> *strIdVector;
    int minlen[2];
    int *isQuery;
// ed
	char dp[STRSIZE][STRSIZE];
	vector<string> strVector;
	vector<pair<string, int> > edShortVector;
	unordered_map<ull, vector<int> *> edMap[STRSIZE][LEN];
    
    bool cmp(const int &x, const int &y);
    int InsertHash(const char *str, int len);
	double ComputeJacc(int lineid, int sz, double threshold);
	void CreateJaccIdf(const char *filename, int id);
	void CreateJaccIndex(const char *filename1, const char *filename2, double threshold);
    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    
	int ComputeEd(const char* str1, int m, const char* str2, int n, int threshold);
    void createEdIndex(const char* filename, int threshold);
    void searchED(const char *querystr, int threshold, vector<pair<int, int>> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold, std::vector<EDJoinResult> &result);
};

#endif
