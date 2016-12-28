#include <string>
#include <algorithm>
//#include <cctype>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <sstream>
#include "strformat.h"

using namespace std;

int main(int argc, char * argv[]) {
	string str = "113a333b 33345 c 981 d 777 o3894)(333.444)";
	vector<int> tmpv = extractInt(str);
	for (int i=0; i<tmpv.size(); i++)
		cout << tmpv[i] << " ";
	cout << endl;

	int xx=1956;
	cout << itos(xx)+",x,xx," << endl;
   //str="abc1234edfG"; this will not work
   str="1234edfG";
   cout << str << " " << stoi(str) << endl;

	return 0;
}
