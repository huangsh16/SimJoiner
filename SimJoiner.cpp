#include "SimJoiner.h"

using namespace std;

SimJoiner::SimJoiner() {
	joinTimes = 0;
	isAppear = new int[300000];
	memset(isAppear, 0, sizeof(isAppear[0]));
	for(int i = 0; i < STRSIZE; ++i) dp[i][0] = i;
    for(int j = 0; j < STRSIZE; ++j) dp[0][j] = j;
    strIdVector = new vector<int> [STRNUM];
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


void SimJoiner::CreateJaccIdf(const char *filename, int id) {
	ifstream fin(filename);
	strSetVector[id].clear();
	char str[STRSIZE]; string strString;
	while(fin.getline(str, STRSIZE)) {
		int last = 0, len = strlen(str); strString = str; 
		set<string> *tmpVec = new set<string>();
		for(int i = 0; i < len; ++i)
            if(str[i] == ' ') {
                if(last < i) 
                    jaccardTrie.Insert(str + last, i - last), tmpVec -> insert(strString.substr(last, i - last));
                last = i + 1;
            }
        if(last < len)
        	jaccardTrie.Insert(str + last, len - last), tmpVec -> insert(strString.substr(last, len - last));
        strSetVector[id].pb(tmpVec);
	}
	fin.close();
}

void SimJoiner::CreateJaccIndex(const char *filename1, const char *filename2, double threshold) {
	// todo : jaccardTrie.clear();
	CreateJaccIdf(filename1, 0);
	CreateJaccIdf(filename2, 1);
	for(int i = 0; i < strSetVector[1].size(); ++i) {
		vector<pair<int, int> > countIdVector;
		for(auto &str : *strSetVector[1][i])
			countIdVector.pb(jaccardTrie.Search(str.c_str(), str.length()));
		sort(countIdVector.begin(), countIdVector.end());
		int preLen = (1 - threshold) * strSetVector[1][i] -> size() + 1;
		for(int j = 0; j < preLen; ++j)
			strIdVector[countIdVector[j].B].pb(i);
	}	
}

double SimJoiner::ComputeJacc(set<string> *l1, set<string> *l2, double threshold) {
    int cnt = 0;
    for (auto& w : *l1) if(l2->find(w) != l2->end()) cnt++;
    return ((double)cnt / (double)(l1->size() + l2->size() - cnt));
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2, double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    CreateJaccIndex(filename1, filename2, threshold);
    for(int i = 0; i < strSetVector[0].size(); ++i) {
    	vector<pair<int, int> > countIdVector;
    	joinTimes++;
		for(auto &str : *strSetVector[0][i])
			countIdVector.pb(jaccardTrie.Search(str.c_str(), str.length()));
		sort(countIdVector.begin(), countIdVector.end());
		int preLen = (1 - threshold) * strSetVector[0][i] -> size() + 1;
		for(int j = 0; j < preLen; ++j)
			for(auto lineId : strIdVector[countIdVector[j].B]) if(isAppear[lineId] != joinTimes){
				isAppear[lineId] = joinTimes;
				double jacc = ComputeJacc(strSetVector[0][i], strSetVector[1][lineId], threshold);
				if(jacc >= threshold)
					result.pb(JaccardJoinResult(i, lineId, jacc));
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
