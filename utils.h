#ifndef __UTIL_H__
#define __UTIL_H__

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
private:
  int cost;
  Direction direction;
  bool isMatch;

public:
  Cell(int value, Direction direction, bool isMatch) : cost(value), direction(direction), isMatch(isMatch) {}
  Cell(const Cell &c) : cost(c.cost), direction(c.direction), isMatch(c.isMatch) {}

  int getCost() const { return cost; }
  void setCost(int value) {cost = value;}
  int getDirection() const {return direction;}
  void setDirection(Direction value) {direction = value;}
  bool getMatch() const {return isMatch;}
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
  int rows, cols;

public:
  Matrix(unsigned rows, unsigned cols) : rows(rows), cols(cols)
  {
    for (unsigned i = 0; i < rows + 1; i++)
    {
      data.push_back(vector<Cell>(cols + 1, Cell(0, NONE, false)));
    }
  }

  Cell &operator()(int row, int col) { return data[row + 1][col + 1]; }
  Cell operator()(int row, int col) const { return data[row + 1][col + 1]; }
  int get_rows() const { return rows; }
  int get_cols() const { return cols; }

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

#endif