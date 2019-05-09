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

struct TrieNode{
	TrieNode* child[128];
	int id, cnt;
	
	TrieNode(): id(0), cnt(0) {
		for(int i = 0; i < 128; ++i)
			child[i] = NULL;
	}
};

struct Trie
{
	TrieNode* root;
	int tot;
	Trie() { 
		root = new TrieNode(); 
		tot = 0;
	}

	void Insert(const char* str, int len) {
		TrieNode* nowNode = root;
		for (int i = 0; i < len; ++i) {
			if(nowNode -> child[(int)str[i]] == NULL) 
				nowNode -> child[(int)str[i]] = new TrieNode();
			nowNode = nowNode -> child[(int)str[i]];
		}
		if(!nowNode -> id) nowNode -> id = ++tot;
		nowNode -> cnt++; 
	}
	pair<int, int> Search(const char* str, int len) {
		TrieNode* nowNode = root;
		for(int i = 0; i < len; ++i) {
			if(nowNode -> child[(int)str[i]] == NULL)
				return mp(-1, -1);
			nowNode = nowNode -> child[(int)str[i]];
		}
		return mp(nowNode -> cnt, nowNode -> id);
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

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();

    int joinTimes;
	int *isAppear;
//jacc
    Trie jaccardTrie;
    vector<set<string>* > strSetVector[2];
    vector<int> strIdVector[STRNUM];
// ed
	int dp[STRSIZE][STRSIZE];
	vector<string> strVector;
	vector<pair<string, int> > edShortVector;
	unordered_map<ull, vector<int> *> edMap[STRSIZE][LEN];
    
	double ComputeJacc(set<string> *l1, set<string> *l2, double threshold);
	void CreateJaccIdf(const char *filename, int id);
	void CreateJaccIndex(const char *filename1, const char *filename2, double threshold);
    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    
	int ComputeEd(const char* str1, int m, const char* str2, int n, int threshold);
    void createEdIndex(const char* filename, int threshold);
    void searchED(const char *querystr, int threshold, vector<pair<int, int>> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold, std::vector<EDJoinResult> &result);
};

#endif
