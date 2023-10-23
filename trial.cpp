#include <bits/stdc++.h>
using namespace std;

int main(){
    vector<int> str{1,2,3,4,5};
    vector<int> str2(str.begin(), str.begin()+3);
    for(auto i: str2){
        cout << i << " ";
    }
}