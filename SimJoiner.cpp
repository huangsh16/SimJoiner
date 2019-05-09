#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
	joinTimes = 0;
	isAppear = new int[300000];
	memset(isAppear, 0, sizeof(isAppear[0]));
	for(int i = 0; i < STRSIZE; ++i) dp[i][0] = i;
    for(int j = 0; j < STRSIZE; ++j) dp[0][j] = j;
}

SimJoiner::~SimJoiner() {
}

#define Abs(_) ((_) > 0 ? (_) : -(_))

unsigned long long HashValue(const char *str, int len) {
    unsigned long long ans = 0;
    for(int i = 0; i < len; ++i)
        ans = ans * 131 + str[i];
    return ans;
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    return SUCCESS;
}

int SimJoiner::ComputeEd(const char* str1, int m, const char* str2, int n, int threshold)
{
    if (Abs(m - n) > threshold)
        return INF;
    for(int i = 1; i <= m; i++) {
        int begin = max(i - threshold, 1);
        int end = min(i + threshold, n);
        if (begin > end)
            break;
        for (int j = begin; j <= end; j++) {
            int t = !(str1[i - 1] == str2[j - 1]);
            int d1 = Abs(i - 1 - j) > threshold ? INF : dp[i - 1][j];
            int d2 = Abs(i - j + 1) > threshold ? INF : dp[i][j - 1];
            dp[i][j] = min(dp[i - 1][j - 1] + t, min(d1 + 1, d2 + 1));
        }
    }
    return dp[m][n];
}

void SimJoiner::createEdIndex(const char* filename, int threshold) {
	strVector.clear();
	edShortVector.clear();
	for(int i = 0; i < STRSIZE; ++i) 
		for(int j = 0; j < LEN; ++j)
			edMap[i][j].clear();

	ifstream fin(filename);
	char str[STRSIZE];
	for(int lineId = 0; fin.getline(str, STRSIZE); ++lineId) {
		int len = strlen(str);
		strVector.pb(str);
		if(len <= threshold) {
			edShortVector.pb(mp(str, lineId));
			continue;
		}
		int len1 = len / (threshold + 1), remain = len - len1 * (threshold + 1);
		int len2 = len1 + (remain > 0), lenPrefix = len - remain * len2;
		for(int i = 0; i < threshold + 1; ++i) {
			ull hashValue;
			if(i * len1 < lenPrefix)hashValue = HashValue(str + i * len1, len1);
			else hashValue = HashValue(str + len - (threshold + 1 - i) * len2, len2);
			auto it = edMap[len][i].find(hashValue);
			if(it == edMap[len][i].end()) {
				vector<int> *tmpVec = new vector<int>();
				tmpVec -> pb(lineId);
				edMap[len][i][hashValue] = tmpVec;
			}
			else it -> B -> pb(lineId);
		}
	}
	fin.close();
}

void SimJoiner::searchED(const char *querystr, int threshold, vector<pair<int, int>> &result){
	joinTimes++;
	result.clear();
	int querylen = strlen(querystr), l = max(0, querylen - threshold), r = min(256, querylen + threshold);
	int len1, remain, len2, lenPrefix, delta, pos, ll, rr, leni; 
	ull hashValue; int ed;
	vector<int> *edMapVector;
	//cout << "searchED querystr : " << querystr << endl;
	for(int len = l; len <= r; ++len) {
		if(edMap[len][0].empty()) continue;
		len1 = len / (threshold + 1), remain = len - len1 * (threshold + 1);
		len2 = len1 + (remain > 0), lenPrefix = len - remain * len2;
		delta = querylen - len;
		for(int i = 0; i < threshold + 1; ++i) {
			if(i * len1 < lenPrefix)pos = i * len1, leni = len1;
			else pos = len - (threshold + 1 - i) * len2, leni = len2;
			ll = max(0, max(pos - i, pos + delta + i - threshold));
			rr = min(querylen - leni, min(pos + i, pos + delta - i + threshold));
			for(int j = ll; j <= rr; ++j) {
				//cout << "find hash [" << j << "," << leni << "]" << endl;
				hashValue = HashValue(querystr + j, leni);
				edMapVector = edMap[len][i][hashValue];
				if(edMapVector) for(auto& lineId : *edMapVector) {
					//cout << "lineId : " << lineId << endl; 
					if(isAppear[lineId] != joinTimes) {
						isAppear[lineId] = joinTimes;
						ed = ComputeEd(querystr, querylen, strVector[lineId].c_str(), strVector[lineId].length(), threshold);
						if(ed <= threshold) result.pb(mp(lineId, ed));
					}
				}
			}	
		}	
	}
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
    result.clear();
    // todo : swap(file1, file2)
    createEdIndex(filename2, threshold);
    ifstream fin(filename1);
    string str;
    std::vector<pair<int, int> > tmpResult;
    for(int i = 0; getline(fin, str); ++i) {
    	searchED(str.c_str(), threshold, tmpResult);
    	for(auto &it : tmpResult)
    		result.pb(EDJoinResult(i, it.A, it.B));
    	for(auto& edShort : edShortVector) {
    		int ed = ComputeEd(str.c_str(), str.length(), edShort.A.c_str(), edShort.A.length(), threshold);
    		if(ed <= threshold)
    			result.pb(EDJoinResult(i, edShort.B, ed));
    	}
    }
    fin.close();
    // todo : sort result with the same i
    sort(result.begin(), result.end());
    return SUCCESS;
}
