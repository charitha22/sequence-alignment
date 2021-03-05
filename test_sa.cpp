#include <iostream>
#include <vector>
#include "needleman_wunsch_sa.h"
#include "smith_waterman_sa.h"
using namespace std;

int match(char a, char b)
{
  if (a == b)
    return 3;
  else
    return -3;
}

int main(int argc, char **argv)
{

  NeedlemanWunschSA<char, vector<char>> NWSA(match);
  SmithWatermanSA<char, vector<char>> SWSA(match);

  string s1 = "TGTTACGG";
  string s2 = "GGTTGACTA";


  vector<char> seq1(s1.begin(), s1.end());
  vector<char> seq2(s2.begin(), s2.end());
  auto solnNW = NWSA.compute(seq1, seq2);
  auto solnSW = SWSA.compute(seq1, seq2);

  cout << "Inputs : \n";
  cout << s1 << "\n" << s2 << "\n";

  cout << "Needleman-Wunsch Solution : \n";
  for (auto pt : solnNW)
    cout << pt << endl;

  cout << "Smith-Waterman solution : \n";
  for (auto pt : solnSW)
    cout << pt << endl;

  return 0;
}