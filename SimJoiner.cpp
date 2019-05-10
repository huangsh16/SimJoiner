#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
	joinTimes = 0;
	isAppear = new int[300000];
	isQuery = new int[38400000];
	memset(isAppear, 0, sizeof(isAppear[0]));
	for(int i = 0; i < STRSIZE; ++i) dp[i][0] = i;
    for(int j = 0; j < STRSIZE; ++j) dp[0][j] = j;
    strIdVector = new vector<int> [STRNUM];
	idToCnt = new int[STRNUM];
	idNum = 0;
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

int SimJoiner::InsertHash(const char *str, int len) {
	ull hashValue = HashValue(str, len);
	int t;
	if(strHashToId[hashValue] == 0) t = strHashToId[hashValue] = ++idNum;
	else t = strHashToId[hashValue];
	idToCnt[t]++;
	return t;
}

void SimJoiner::CreateJaccIdf(const char *filename, int id) {
	ifstream fin(filename);
	idVectorVector[id].clear();
	minlen[id] = STRSIZE;
	char str[STRSIZE];
	while(fin.getline(str, STRSIZE)) {
		int last = 0, len = strlen(str); 
		vector<int> tmpVec;
		for(int i = 0; i < len; ++i)
            if(str[i] == ' ') {
                if(last < i) 
                    tmpVec.pb(InsertHash(str + last, i - last));
                last = i + 1;
            }
        if(last < len)
			tmpVec.pb(InsertHash(str + last, len - last));
		sort(tmpVec.begin(), tmpVec.end());
		tmpVec.resize(unique(tmpVec.begin(), tmpVec.end()) - tmpVec.begin());
		minlen[id] = min(minlen[id], (int)tmpVec.size());
        idVectorVector[id].pb(tmpVec);
	}
	fin.close();
}

void SimJoiner::CreateJaccIndex(const char *filename1, const char *filename2, double threshold) {
	// todo : jaccardTrie.clear();
	CreateJaccIdf(filename1, 0);
	CreateJaccIdf(filename2, 1);
	for(int i = 0; i < idVectorVector[1].size(); ++i) {
		vector<pair<int, int> > tmpVec;
		for(int j = 0; j < idVectorVector[1][i].size(); ++j)
			tmpVec.pb(mp(idToCnt[idVectorVector[1][i][j]], idVectorVector[1][i][j]));
		sort(tmpVec.begin(), tmpVec.end());
		TmpVec.pb(tmpVec);
		int preLen = 1 + 1.0 / (1 + threshold) * idVectorVector[1][i].size() - threshold / (1 + threshold) * minlen[0];
		for(int j = 0; j < preLen; ++j)
			strIdVector[tmpVec[j].B].pb(i);
	}	
}

double SimJoiner::ComputeJacc(int lineId, int sz, double threshold) {
	int cnt = 0, sz1 = TmpVec[lineId].size();
	for(int i = 0; i < sz1; ++i) {
		if(isQuery[TmpVec[lineId][i].B] == joinTimes) ++cnt;
		if((cnt + sz1 - i - 1) < threshold * (sz - cnt + i + 1)) return 0;
	}
	return 1.0 * cnt / (sz + sz1 - cnt);
}

pair<int, int> Count(vector<pair<int, int> > &v1, int x, vector<pair<int, int> > &v2, int y) {
	int i = 0, j = 0, cnt = 0;
	for(; i <= x; ++i) {
		while(v2[j] < v1[i]) ++j;
		cnt += v2[j] == v1[i];
	}
	return mp(cnt, x + y + 2 - cnt);
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    CreateJaccIndex(filename1, filename2, threshold);
    for(int i = 0; i < idVectorVector[0].size(); ++i) {
    	joinTimes++;
    	for(auto id : idVectorVector[0][i]) isQuery[id] = joinTimes;
		int preLen = 1 + 1.0 / (1 + threshold) * idVectorVector[0][i].size() - threshold / (1 + threshold) * minlen[1];

		vector<pair<int, int> > tmpVec;
		for(int j = 0; j < idVectorVector[0][i].size(); ++j)
			tmpVec.pb(mp(idToCnt[idVectorVector[0][i][j]], idVectorVector[0][i][j]));
		sort(tmpVec.begin(), tmpVec.end());
		for(int j = 0; j < preLen; ++j) for(auto lineId : strIdVector[tmpVec[j].B]) {
			if(isAppear[lineId] != joinTimes){
				isAppear[lineId] = joinTimes;
/*
				int px = 1 + 1.0 / (1 + threshold) * idVectorVector[0][i].size() - threshold / (1 + threshold) * idVectorVector[1][lineId].size();
				int py = 1 + 1.0 / (1 + threshold) * idVectorVector[1][lineId].size() - threshold / (1 + threshold) * idVectorVector[0][i].size();
				if(px && py) {
					int x = px - 1, y = py - 1;
					while(tmpVec[x] != TmpVec[lineId][y]) {
						while(tmpVec[x] > TmpVec[lineId][y]) x--;
						while(tmpVec[x] < TmpVec[lineId][y]) y--;
					}
					auto tt = Count(tmpVec, x, TmpVec[lineId], y);
					int ubound = tt.A + min(tmpVec.size() - x - 1, TmpVec[lineId].size() - y - 1);
					int lbound = tt.B + max(tmpVec.size() - x - 1, TmpVec[lineId].size() - y - 1);
					if(ubound < threshold * lbound) continue;
				}
*/
				double jacc = ComputeJacc(lineId, idVectorVector[0][i].size(), threshold);
				if(jacc >= threshold)
					result.pb(JaccardJoinResult(i, lineId, jacc));
			}
		}
	}
    sort(result.begin(), result.end());
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
				hashValue = HashValue(querystr + j, leni);
				edMapVector = edMap[len][i][hashValue];
				if(edMapVector) for(auto& lineId : *edMapVector) if(isAppear[lineId] != joinTimes) {
					isAppear[lineId] = joinTimes;
					ed = ComputeEd(querystr, querylen, strVector[lineId].c_str(), strVector[lineId].length(), threshold);
					if(ed <= threshold) result.pb(mp(lineId, ed));
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
    // done : sort result with the same i
    sort(result.begin(), result.end());
    return SUCCESS;
}
