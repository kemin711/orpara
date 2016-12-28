#include <string>
#include <algorithm>
//#include <cctype>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <sstream>
#include "strformat.h"
#include "stddev.h"

using namespace std;

//using namespace KZUtility;
void teststddev() {
   cout << "testing stddev lib\n";
   stddev stat;
   double values[11] = { 0.5, 2.4, 3.7, 4.2, 2.8, 3.1, 2.1, 1.9, 3.4, 4.1, 2.9 };
   for (int i=0; i<11; ++i) {
      stat(values[i]);
   }
   cout << "Mean: " << stat.getMean() << "  stddev: "
      << stat.getStd() << endl;
   string x20(35, 'x');
   cout << "done with stddev test\n" << x20 << "\n";
   
}
void teststddevN() {
   cout << "testing stddev lib with multiple values\n";
   stddev stat;
   double values[11] = { 0.5, 2.4, 3.7, 4.2, 2.8, 3.1, 2.1, 1.9, 3.4, 4.1, 2.9 };
   int count[11] = {1, 3, 7, 2, 3, 4, 8, 7, 9, 4, 5};
   for (int i=0; i<11; ++i) {
      cout << values[i] << '(' << count[i] << ')';
      if (i < 10) cout << ", ";
      stat(values[i], count[i]);
   }
   cout << "\nMean: " << stat.getMean() << "  stddev: "
      << stat.getStd() << endl;
   cout << "done with stddev test\n" << string(40, 'y') << "\n";
   
}

int main(int argc, char * argv[]) {
   teststddev();
   teststddevN();
	string str = "113a333b 33345 c 981 d 777 o3894)(333.444)";
	cout << "total words of " << str << " is " << wc(str) << endl;
	vector<int> tmpv = extractInt(str);
	for (int i=0; i<tmpv.size(); i++)
		cout << tmpv[i] << " ";
	cout << endl;

	int xx=1956;
	cout << itos(xx)+",x,xx," << endl;
	str = "\tAAA\t345\tBBBBB\t\t";
	cout << "trying to split str\n" << str << endl;
	vector<string> tmpvs = split(str, '\t');
	for (xx=0; xx<tmpvs.size(); xx++) {
		cout << "element " << xx << " is " << tmpvs[xx] << endl;
	}
	str = "long phase with several words";
	cout << str << "\nhas acronym:\n" << acronym(str)
		<< " and acronym(2)" << acronym(str, 2) << endl;
	str = "16pHQG;4 24b2/STAC2 APOC1'";
	string str1 = delall(str, ";/'");
	cout << str << "\nafter cleaning\n" << str1 << endl;
	str = "with      multiple \t spaces    lslsl\n";
	cout << "str before singleSpace()\n" << str << endl;
	singleSpace(str);
	cout << str << endl;
	str = "NPY Y4";
	cout << "acronymWithTag(): " << str << " => " << acronymWithTag(str) 
		<< endl;
	str = "1.3,4-5";
	if (isnumber(str, ".,-")) {
		cout << str << " is number\n";
	}
	else cout << str << " is not number\n";

	return 0;
}


