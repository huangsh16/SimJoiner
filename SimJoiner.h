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


struct TrieNode{
	TrieNode* child[128];
	vector<int> indexVector;
	TrieNode() {
		for(int i = 0; i < 128; ++i)
			child[i] = NULL;
	}
};

struct Trie
{
	TrieNode* root;
	Trie() { 
		root = new TrieNode(); 
	}
	void Insert(const char* str, int len, int lineId) {
		TrieNode* nowNode = root;
		for (int i = 0; i < len; ++i) {
			if(nowNode -> child[(int)str[i]] == NULL) 
				nowNode -> child[(int)str[i]] = new TrieNode();
			nowNode = nowNode -> child[(int)str[i]];
		}
		if(nowNode -> indexVector.empty() || *(nowNode -> indexVector.end() - 1) != lineId)
			nowNode -> indexVector.push_back(lineId);
	}
	vector<int>* Search(const char* str, int len) {
		TrieNode* nowNode = root;
		for(int i = 0; i < len; ++i) {
			if(nowNode -> child[(int)str[i]] == NULL)
				return NULL;
			nowNode = nowNode -> child[(int)str[i]];
		}
		return &(nowNode -> indexVector);
	}
};

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

const int SUCCESS = 0;
const int FAILURE = 1;
const int STRSIZE = 260;
const int LEN = 10;
const int INF = (int)1e9;
#define pb push_back
#define A first
#define B second
#define mp make_pair
typedef unsigned long long ull;

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();

    int joinTimes;
// ed
	int *isAppear;
	int dp[STRSIZE][STRSIZE];
	vector<string> strVector;
	vector<pair<string, int> > edShortVector;
	unordered_map<ull, vector<int> *> edMap[STRSIZE][LEN];
    

    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    
	int ComputeEd(const char* str1, int m, const char* str2, int n, int threshold);
    void createEdIndex(const char* filename, int threshold);
    void searchED(const char *querystr, int threshold, vector<pair<int, int>> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold, std::vector<EDJoinResult> &result);
};

#endif
