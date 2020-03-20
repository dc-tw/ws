#include <iostream>
#include <string>
//using namespace std;

string func(std::string str, int len, bool strand){
    std::string ans;
    if(strand){
        for(int i=0; i<len; ++i){
            switch(str[i]){
                case 'A':
                    ans = ans + 'A';
                break;
                case 'T':
                    ans = ans + 'T';
                break;
                case 'C':
                    ans = ans + 'C';
                break;
                case 'G':
                    ans = ans + 'G';
                break;
                case 'Z':
                    ans = ans + 'T';
                break;
                case 'Q':
                    ans = ans + 'G';
                break;
                default:
                break;
            }
        }
    }
    else{
        for(int i=0; i<len; ++i){
            switch(str[i]){
                case 'A':
                    ans = ans + 'A';
                break;
                case 'T':
                    ans = ans + 'T';
                break;
                case 'C':
                    ans = ans + 'C';
                break;
                case 'G':
                    ans = ans + 'G';
                break;
                case 'Z':
                    ans = ans + 'C';
                break;
                case 'Q':
                    ans = ans + 'A';
                break;
                default:
                break;
            }
        }
    }
    return ans;
}

string func2(std::string str, int len, bool strand){
    std::string ans;
    char arr1[26], arr2[26];
    arr1['A'-'A']='A';arr1['T'-'A']='T';arr1['C'-'A']='C';arr1['G'-'A']='G';arr1['Z'-'A']='T';arr1['Q'-'A']='G';
    arr2['A'-'A']='A';arr2['T'-'A']='T';arr2['C'-'A']='C';arr2['G'-'A']='G';arr2['Z'-'A']='C';arr2['Q'-'A']='A';
    if(strand){
        for(int i=0; i<len; ++i){
            ans += arr1[str[i]-'A'];
        }
    }
    else{
        for(int i=0; i<len; ++i){
            ans += arr2[str[i]-'A'];
        }
    }
    return ans;
}

int main() {
    std::string str = "ACGZTQA", ans;
    int len=7;
    ans = func(str, len, 0);
   
    cout<<ans;
}