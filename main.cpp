#include "SimJoiner.h"

using namespace std;

int main(int argc, char **argv) {
    SimJoiner joiner;

    vector<EDJoinResult> resultED;
    vector<JaccardJoinResult> resultJaccard;

    unsigned edThreshold = 2;
    double jaccardThreshold = 0.85;

    cout << "yes\n";

//    joiner.joinJaccard(argv[1], argv[2], jaccardThreshold, resultJaccard);
//    joiner.joinED(argv[1], argv[2], edThreshold, resultED);

    cout << "empty :" << resultJaccard.empty() << endl;

    for(auto i : resultJaccard) 
        cout << i.id1 << " " << i.id2 << " " << i.s << endl;


/*  cout << "empty :" << resultED.empty() << endl;

    for(auto i : resultED) 
    	cout << i.id1 << " " << i.id2 << " " << i.s << endl;

    joiner.joinED(argv[2], argv[1], edThreshold, resultED);

    cout << "empty :" << resultED.empty() << endl;

    for(auto i : resultED) 
        cout << i.id1 << " " << i.id2 << " " << i.s << endl;
*/
    return 0;
}
