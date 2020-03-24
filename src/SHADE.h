/* 
  SHADE.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _SHADE_H_
#define _SHADE_H_

#include "DE.h"

class SHADE: public DE {
 public:
  using DE::DE;
  virtual void run();
};

#endif
