#ifndef RANDOMBASE_H
#define RANDOMBASE_H

#include <random>

using namespace std;

namespace orpara {

/**
 * Ambiguous base 
 * code         Complement
 * A                  T
 * C                  G
 * G                  C
 * T                  A
 * Y: Pyrimidine [CT] R
 * R: Purine [AG]     Y
 * W: weak [AT]       W
 * S: strong [CG]     S
 * K: keto [GT]       M         
 * M: amino [AC]      K        
 * D: [AGT]           H 
 * V: [ACG]           B 
 * H: [ACT]           D 
 * B: [CGT]           V 
 * X/N  any           X/N
 * -  Gap             -
 */
class RandomBase {
   public:
      static RandomBase& getInstance() {
         static RandomBase single;
         return single;
      }
      RandomBase(const RandomBase& rb) = delete;
      RandomBase& operator=(const RandomBase& rb) = delete;
      ~RandomBase() { }
      /**
       * @return a random base from A,C,G, or T
       */
      char operator()() {
         return baseN[udis4(mteng)];
      }
      char operator()(char ambiguous) {
         char AM=toupper(ambiguous);
         switch(AM) {
            case 'N' : case 'X': return baseN[udis4(mteng)];
            case 'Y': return baseY[udis2(mteng)];
            case 'R': return baseR[udis2(mteng)];
            case 'W': return baseW[udis2(mteng)];
            case 'S': return baseS[udis2(mteng)];
            case 'K': return baseK[udis2(mteng)]; 
            case 'M': return baseM[udis2(mteng)];
            case 'D': return baseD[udis3(mteng)];
            case 'V': return baseV[udis3(mteng)];
            case 'H': return baseH[udis3(mteng)];
            case 'B': return baseB[udis3(mteng)];
            default: throw runtime_error("wrong base: " + string(1, AM));
         }
      }

   private:
      RandomBase() 
         : randev(), mteng(randev()), udis4(0,3), udis3(0,2), udis2(0,1),
           baseN{'A', 'C', 'G', 'T'},
           baseY{'C', 'T'}, baseR{'A', 'G'},
           baseW{'A', 'T'}, baseS{'C', 'G'},
           baseK{'G', 'T'}, baseM{'A', 'C'},
           baseD{'A', 'G', 'T'}, baseV{'A', 'C', 'G'},
           baseH{'A', 'C', 'T'}, baseB{'C', 'G', 'T'}
      { }
      std::random_device randev;
      std::mt19937_64 mteng;
      std::uniform_int_distribution<> udis4;
      std::uniform_int_distribution<> udis3;
      std::uniform_int_distribution<> udis2;
      char baseN[4];
      char baseY[2];
      char baseR[2];
      char baseW[2];
      char baseS[2];
      char baseK[2];
      char baseM[2];
      char baseD[3];
      char baseV[3];
      char baseH[3];
      char baseB[3];
};
}
#endif
