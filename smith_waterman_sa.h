#ifndef __SMITH_WATERMAN_H__
#define __SMITH_WATERMAN_H__
#include <functional>
#include <memory>
#include <ostream>
#include <assert.h>
#include <string>
#include "utils.h"
using namespace std;


template <typename elemTy, typename ArrayTy>
class SmithWatermanSA
{
private:
  function<int(elemTy, elemTy)> match;
  const unsigned gap {2};

  AlignedSeq<elemTy> constructSoln(Matrix &m, ArrayTy &seq1, ArrayTy &seq2, int startRow, int endRow, int startCol, int endCol)
  {
    // find the max cost
    int maxi = endRow, maxj = endCol;
    unsigned maxCost = m(endRow, endCol).getCost();
    for (int i = startRow; i <= endRow; ++i)
    {
      for (int j = startCol; j <= endCol; ++j)
      {
        if (m(i, j).getCost() > maxCost)
        {
          maxi = i;
          maxj = j;
          maxCost = m(i, j).getCost();
        }
      }
    }
    AlignedSeq<elemTy> soln;
    int starti = maxi;
    int startj = maxj;
    while (starti >= startRow && startj >= startCol)
    {
      if (m(starti, startj).getDirection() == DIAG)
      {
        AlignedPair<elemTy> point(&seq1[starti], &seq2[startj]);
        soln.insert(soln.begin(), point);
        starti--;
        startj--;
      }
      else if (m(starti, startj).getDirection() == LEFT)
      {
        AlignedPair<elemTy> point(nullptr, &seq2[startj]);
        soln.insert(soln.begin(), point);
        startj--;
      }
      else
      {
        AlignedPair<elemTy> point(&seq1[starti], nullptr);
        soln.insert(soln.begin(), point);
        starti--;
      }
    }

    // find the reminder solns from end of the Matrix
    if (maxi < endRow && maxj < endCol)
    {
      auto solnRight = constructSoln(m, seq1, seq2, maxi + 1, endRow, maxj + 1, endCol);
      soln.concat(solnRight);
    }

    // add any remainder from the begining of the Matrix
    while (starti >= startRow)
    {
      AlignedPair<elemTy> point(&seq1[starti], nullptr);
      soln.insert(soln.begin(), point);
      starti--;
    }
    while (startj >= startCol)
    {
      AlignedPair<elemTy> point(nullptr, &seq2[startj]);
      soln.insert(soln.begin(), point);
      startj--;
    }

    return soln;
  }

public:
  SmithWatermanSA(function<int(elemTy, elemTy)> match) : match(match) {}

  AlignedSeq<elemTy> compute(ArrayTy &seq1, ArrayTy &seq2)
  {
    Matrix m(seq1.size(), seq2.size());

    for (unsigned i = 0; i < m.get_rows(); ++i)
    {
      for (unsigned j = 0; j < m.get_cols(); ++j)
      {
        int diag_cost = m(i - 1, j - 1).getCost() + match(seq1[i], seq2[j]);
        int left_cost = m(i, j - 1).getCost() - gap;
        int top_cost = m(i - 1, j).getCost() - gap;
        int cost = diag_cost;
        Direction d = DIAG;
        bool isMatch = true;

        if (cost < left_cost)
        {
          d = LEFT;
          cost = left_cost;
          isMatch = false;
        }

        if (cost < top_cost)
        {
          d = TOP;
          cost = top_cost;
          isMatch = false;
        }

        if (cost < 0)
        {
          cost = 0;
        }

        m(i, j) = Cell(cost, d, isMatch);
      }
    }

    return constructSoln(m, seq1, seq2, 0, m.get_rows() - 1, 0, m.get_cols() - 1);
  }
};

#endif