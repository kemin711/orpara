#include <iostream>
#include <assert.h>

#include <range.h>

int main(int argc, char* argv[]) {
   Range rg;
   cout << "Default range: " << rg << endl;
   Range largerg(-100, 999999);
   cout << "Large range: " << largerg << endl;
   assert(rg.length());
   Range rg2(3, 87);
   Range rg3(8, 200);

   if (rg2 < rg3) {
      cout << rg2 << " is smaller than " << rg3 << endl;
   }
   else {
      cout << rg2 << " is not smaller than " << rg3 << endl;
   }

   return 0;
}
