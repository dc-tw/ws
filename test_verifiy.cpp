#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <cstdlib>

using namespace std;

void function(string str, int len, int prob, int quality) {
	char arr1[26], arr2[26];
    arr1['A'-'A']='A';arr1['T'-'A']='T';arr1['C'-'A']='C';arr1['G'-'A']='G';arr1['Z'-'A']='T';arr1['Q'-'A']='G';
    arr2['A'-'A']='A';arr2['T'-'A']='T';arr2['C'-'A']='C';arr2['G'-'A']='G';arr2['Z'-'A']='C';arr2['Q'-'A']='A';
	
    /*---determin length---*/
    int start = rand() % len;
    while(start>10)start = rand() % len;
    int end = rand() % len;
    if(end < start)end = len;

    /*---scan string---*/
    string ans;
    int strand = rand()%2;
    for(int i=start; i<=end; ++i){
        if(strand){
	        ans += arr1[str[i]-'A'];
	    }
	    else{
	        ans += arr2[str[i]-'A'];
	    }
    }

    cout<<"start : "<< start<<" end : "<< end<<endl;
    cout<<ans<<endl;
}



int main() {
    srand(time(NULL));

    string s = "AATCGTACTCAAGTT";
    int tmp, prob = 5;
    for(int i=0; i<s.length(); ++i){
        if(s[i]=='C'){
            tmp = rand()%10;
            if(prob<=tmp)s[i] = 'Q';
        }
        else if(s[i]=='G'){
            tmp = rand()%10;
            if(prob<=tmp)s[i] = 'Z';
        }
    }
    cout<<s<<endl;
    for(int j=0; j<3; ++j)
		function(s, s.length(), 5, 5);
    
    return 0;
}
