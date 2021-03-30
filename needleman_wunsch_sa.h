#ifndef __NEEDLEMAN_WUNSCH_H__
#define __NEEDLEMAN_WUNSCH_H__
#include <functional>
#include <memory>
#include <ostream>
#include <assert.h>
#include <string>
#include "utils.h"

using namespace std;

template <typename elemTy, typename ArrayTy>
class NeedlemanWunschSA
{
private:
  ScoringFunction<elemTy>& costModel;

  AlignedSeq<elemTy> constructSoln(Matrix &m, ArrayTy &seq1, ArrayTy &seq2)
  {
    AlignedSeq<elemTy> result;

    int i = m.get_rows() - 1;
    int j = m.get_cols() - 1;

    while (i >= 0 || j >= 0)
    {

      if (i >= 0 && j >= 0 && m(i, j).getMatch())
      {
        AlignedPair<elemTy> point(&seq1[i], &seq2[j]);
        result.insert(result.begin(), point);
        i--;
        j--;
      }
      else if (i >= 0 && m(i, j).getDirection() == TOP)
      {
        AlignedPair<elemTy> point(&seq1[i], nullptr);
        result.insert(result.begin(), point);
        i--;
      }
      else
      {
        AlignedPair<elemTy> point(nullptr, &seq2[j]);
        result.insert(result.begin(), point);
        j--;
      }
    }

    return result;
  }

public:
  NeedlemanWunschSA(ScoringFunction<elemTy>& costModel) : costModel(costModel) {}

  AlignedSeq<elemTy> compute(ArrayTy &seq1, ArrayTy &seq2)
  {
    Matrix m(seq1.size(), seq2.size());

    for (int i = -1; i < m.get_rows(); i++)
    {
      m(i, -1).setCost(-costModel.gap(0));
      m(i, -1).setDirection(TOP);
    }
    for (int i = -1; i < m.get_cols(); i++)
    {
      m(-1, i).setCost(-costModel.gap(0));
      m(-1, i).setDirection(LEFT);
    }

    for (unsigned i = 0; i < m.get_rows(); ++i)
    {
      for (unsigned j = 0; j < m.get_cols(); ++j)
      {
        int diag_cost = m(i - 1, j - 1).getCost() + costModel(*(seq1.begin() + i), *(seq2.begin() + j));
        int left_cost = m(i, j - 1).getCost() - costModel.gap(0);
        int top_cost = m(i - 1, j).getCost() - costModel.gap(0);
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

        m(i, j) = Cell(cost, d, isMatch);
      }
    }

    cout << m << "\n";

    return constructSoln(m, seq1, seq2);
  }
};

#endif