/* $Id */
#ifndef _TEST_CLASS_H_
#define _TEST_CLASS_H_

#include <iostream>

class testClass {

public:

	testClass(int a);
	void Print(ostream& out);

private:

	int a_;

};

#endif /* _TEST_CLASS_H_ */