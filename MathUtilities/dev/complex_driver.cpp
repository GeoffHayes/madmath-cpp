/*
 * Development test driver for the TComplexNumber class.
 *
 * 2013-06-01   Geoff Hayes     Initial Release.
 */

#include "TComplexNumber.h"

#if TCOMPLEXNUMBER_TEST_DRIVER

int main()
{
    TComplexNumber<float> a(2,3);
    TComplexNumber<float> b(1,-6);
    TComplexNumber<float> c = a + b;
    c.write();

    a.setReal(5);
    a.setImag(-2);
    b.setReal(-4);
    b.setImag(-1);
    c = a - b;
    c.write();

    a.setReal(2);
    a.setImag(-1);
    b.setReal(3);
    b.setImag(4);
    c = a*b;
    c.write();

    a.setReal(3);
    a.setImag(0);
    b.setReal(0);
    b.setImag(2);
    c = a/b;
    c. write();

    b.setReal(2);
    b.setImag(1);
    c = a/b;
    c.write();

    return 0;
}

#endif
