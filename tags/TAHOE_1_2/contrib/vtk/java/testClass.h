/* $Id */
#ifndef _TEST_CLASS_H_
#define _TEST_CLASS_H_

#include <iostream.h>
#include "iConsoleT.h"

class testClass {

public:

	testClass(int a);
	void Print(ostream& out);
	void SetA (int x);
	int GetA();

	

private:

	int a_;

};

#endif /* _TEST_CLASS_H_ */
