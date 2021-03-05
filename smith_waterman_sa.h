#ifndef __SMITH_WATERMAN_H_
#define __SMITH_WATERMAN_H_
#include <functional>
#include <memory>
#include <ostream>
#include <assert.h>
#include <string>
using namespace std;

enum Direction
{
  LEFT,
  TOP,
  DIAG,
  NONE
};

struct Cell
{
  Cell(int value, Direction direction, bool isMatch) : cost(value), direction(direction), isMatch(isMatch) {}
  int cost;
  Direction direction;
  bool isMatch;
  friend ostream &operator<<(ostream &os, const Cell &c)
  {
    os << "(" << c.cost << ",";
    if (c.direction == DIAG)
      os << "DIAG";
    else if (c.direction == TOP)
      os << "TOP";
    else if (c.direction == LEFT)
      os << "LEFT";
    else
      os << "NONE";
    os << ",";
    if (c.isMatch)
      os << "T";
    else
      os << "F";
    os << ")";

    return os;
  }
};

class Matrix
{
private:
  vector<vector<Cell>> data;
  unsigned rows, cols;

public:
  Matrix(unsigned rows, unsigned cols) : rows(rows), cols(cols)
  {
    for (unsigned i = 0; i < rows + 1; i++)
    {
      data.push_back(vector<Cell>(cols + 1, Cell(0, NONE, false)));
    }
  }

  Cell &operator()(unsigned row, unsigned col) { return data[row + 1][col + 1]; }
  Cell operator()(unsigned row, unsigned col) const { return data[row + 1][col + 1]; }
  unsigned get_rows() const { return rows; }
  unsigned get_cols() const { return cols; }

  friend ostream &operator<<(ostream &os, const Matrix &m)
  {
    for (unsigned i = 0; i < m.get_rows(); ++i)
    {
      for (unsigned j = 0; j < m.get_cols(); ++j)
      {
        os << m(i, j);
      }
      os << "\n";
    }
    return os;
  }
};

template <typename elemTy>
struct AlignedPair
{
  elemTy *eL = nullptr;
  elemTy *eR = nullptr;
  bool match = false;

public:
  AlignedPair(elemTy *eL, elemTy *eR) : eL(eL), eR(eR), match(eL && eR && *eL == *eR) {}
  bool isMatch() { return match; }
  bool isMisMatch() { return !match; }
  elemTy *getLeft() const { return eL; }
  elemTy *getRight() const { return eR; }
  friend ostream &operator<<(ostream &os, const AlignedPair<elemTy> &ap)
  {
    if (ap.getLeft() != nullptr)
      os << *ap.getLeft();
    else
      os << "_";
    os << " : ";
    if (ap.getRight() != nullptr)
      os << *ap.getRight();
    else
      os << "_";
    return os;
  }
};

template <typename elemTy>
class AlignedSeq : public vector<AlignedPair<elemTy>>
{

public:
  using vector<AlignedPair<elemTy>>::size;
  using vector<AlignedPair<elemTy>>::end;
  using vector<AlignedPair<elemTy>>::begin;
  using vector<AlignedPair<elemTy>>::reserve;
  using vector<AlignedPair<elemTy>>::insert;

  AlignedSeq &concat(AlignedSeq &other)
  {

    reserve(size() + other.size());
    insert(end(), other.begin(), other.end());
    return *this;
  }
};

template <typename elemTy, typename ArrayTy>
class SmithWatermanSA
{
private:
  function<int(elemTy, elemTy)> match;
  unsigned gap = 2;
  int count = 0;
  AlignedSeq<elemTy> constructSoln(Matrix &m, ArrayTy &seq1, ArrayTy &seq2, int startRow, int endRow, int startCol, int endCol)
  {
    // find the max cost
    int maxi = endRow, maxj = endCol;
    unsigned maxCost = m(endRow, endCol).cost;
    for (int i = startRow; i <= endRow; ++i)
    {
      for (int j = startCol; j <= endCol; ++j)
      {
        if (m(i, j).cost > maxCost)
        {
          maxi = i;
          maxj = j;
          maxCost = m(i, j).cost;
        }
      }
    }
    AlignedSeq<elemTy> soln;
    int starti = maxi;
    int startj = maxj;
    while (starti >= startRow && startj >= startCol)
    {
      if (m(starti, startj).direction == DIAG)
      {
        AlignedPair<elemTy> point(&seq1[starti], &seq2[startj]);
        soln.insert(soln.begin(), point);
        starti--;
        startj--;
      }
      else if (m(starti, startj).direction == LEFT)
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

  vector<AlignedPair<elemTy>> compute(ArrayTy &seq1, ArrayTy &seq2)
  {
    Matrix m(seq1.size(), seq2.size());

    for (unsigned i = 0; i < m.get_rows(); ++i)
    {
      for (unsigned j = 0; j < m.get_cols(); ++j)
      {
        int diag_cost = m(i - 1, j - 1).cost + match(seq1[i], seq2[j]);
        int left_cost = m(i, j - 1).cost - gap;
        int top_cost = m(i - 1, j).cost - gap;
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
    cout << m;

    return constructSoln(m, seq1, seq2, 0, m.get_rows() - 1, 0, m.get_cols() - 1);
  }
};

#endif