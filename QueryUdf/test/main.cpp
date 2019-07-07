/******************************************************************************
 * Copyright (c)  2016, TigerGraph Inc.
 * All rights reserved
 ******************************************************************************/
#include <assert.h>
#include "../ExprFunctions.hpp"

/**
 *  Unit testing of the expression functions
 */ 
int main(){
  assert(UDIMPL::str_to_int("12345") == 12345);
  std::cout << "Passed Function str_to_int" << std::endl;

  assert(UDIMPL::float_to_int(123.45) == 123);
  std::cout << "Passed Function float_to_int" << std::endl;

  assert(UDIMPL::to_string(123.45) == "123.45");
  std::cout << "Passed Function to_string" << std::endl;

  return 0;
}

