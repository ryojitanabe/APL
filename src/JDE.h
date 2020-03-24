/* 
  JDE.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _JDE_H_
#define _JDE_H_

#include "DE.h"

class JDE: public DE {
 public:
  using DE::DE;
  virtual void run();
};

#endif
