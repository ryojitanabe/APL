/* 
  ClassicalDE.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _CLASSICALDE_H_
#define _CLASSICALDE_H_

#include "DE.h"

class classicalDE: public DE {
 public:
  using DE::DE;
  virtual void run();
};

#endif
