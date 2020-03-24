/* 
  JADE.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _JADE_H_
#define _JADE_H_

#include "DE.h"

class JADE: public DE {
 public:
  using DE::DE;
  virtual void run();
};

#endif
