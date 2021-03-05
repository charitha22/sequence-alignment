#include <iostream>
#include <vector>
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

  SmithWatermanSA<char, vector<char>> SWSA(match);

  // example from https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
  vector<char> seq2 = {'G', 'G', 'T', 'T', 'G', 'A', 'C', 'T', 'A'};
  vector<char> seq1 = {'T', 'G', 'T', 'T', 'A', 'C', 'G', 'G'};

  auto soln = SWSA.compute(seq1, seq2);

  for (auto &pt : soln)
    cout << pt << "\n";

  return 0;
}