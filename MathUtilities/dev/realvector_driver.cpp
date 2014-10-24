/*
 * Development test driver for the TRealVector class.
 *
 * 2013-06-01   Geoff Hayes     Initial Release.
 */

#include "TRealVector.h"
#include <cstring>
#include <vector>
#include <algorithm>
#include <ostream>

#if TREALVECTOR_TEST_DRIVER

int main()
{
    unsigned int numElements = 10;

    float a[numElements];
    float b[numElements];
    float c[numElements];
    float d[numElements];

    memset(a, 0, numElements*sizeof(float));
    memset(b, 0, numElements*sizeof(float));
    memset(c, 0, numElements*sizeof(float));
    memset(d, 0, numElements*sizeof(float));

    TRealVector<float, unsigned int>::vset(numElements, a, 1);
    TRealVector<float, unsigned int>::vset(numElements, b, 2);
    TRealVector<float, unsigned int>::vset(numElements, c, 3);
    TRealVector<float, unsigned int>::vset(numElements, d, 4);

    TRealVector<float, unsigned int>::vwrite(numElements, a);
    TRealVector<float, unsigned int>::vwrite(numElements, b);
    TRealVector<float, unsigned int>::vwrite(numElements, c);
    TRealVector<float, unsigned int>::vwrite(numElements, d);

    TRealVector<float, unsigned int>::vadd(numElements, a, b);
    TRealVector<float, unsigned int>::vwrite(numElements, a);

    TRealVector<float, unsigned int>::vsub(numElements, a, b);
    TRealVector<float, unsigned int>::vwrite(numElements, a);

    TRealVector<float, unsigned int>::vmult(numElements, a, c);
    TRealVector<float, unsigned int>::vwrite(numElements, a);

    TRealVector<float, unsigned int>::vdiv(numElements, a, c);
    TRealVector<float, unsigned int>::vwrite(numElements, a);

    const float dotprod = TRealVector<float, unsigned int>::vdot(numElements,a,a);
    std::cout << "dotproduct of a*a=" << dotprod << std::endl;

    TRealVector<float, unsigned int>::vlincomb(numElements,d,b,2);
    TRealVector<float, unsigned int>::vwrite(numElements, d);

    TRealVector<float, unsigned int>::vlincomb(numElements,a,a,2);
    TRealVector<float, unsigned int>::vwrite(numElements, a);

    TRealVector<float, unsigned int>::vlincomb(numElements,d,a,b,2,3);
    TRealVector<float, unsigned int>::vwrite(numElements, d);

    TRealVector<float, unsigned int>::vlincomb(numElements,d,a,b,c, 2,3,4);
    TRealVector<float, unsigned int>::vwrite(numElements, d);

    TRealVector<float, unsigned int>::vset(numElements, a, 1);
    const float sum = TRealVector<float, unsigned int>::vsum(numElements, a);
    std::cout << "sum of a=" << sum << std::endl;

    TRealVector<float, unsigned int>::vset(numElements, a, 2);
    const float sumsq = TRealVector<float, unsigned int>::vsumsq(numElements, a);
    std::cout << "sum sq of a=" << sumsq << std::endl;

    const float mean = TRealVector<float, unsigned int>::vmean(numElements, a);
    std::cout << "mean of a=" << mean << std::endl;

    float e[] = {1,5,9};
    float f[] = {7,15,22};

    float stddev = TRealVector<float, unsigned int>::vstddev(3, e);
    std::cout << "stddev of e=" << stddev << std::endl;

    stddev = TRealVector<float, unsigned int>::vstddev(3, f);
    std::cout << "stddev of f=" << stddev << std::endl;

    float g[] = {4, 9};
    float h[] = {-2, 5};

    stddev = TRealVector<float, unsigned int>::vstddev(2, g);
    std::cout << "stddev of g=" << stddev << std::endl;

    stddev = TRealVector<float, unsigned int>::vstddev(2, h);
    std::cout << "stddev of h=" << stddev << std::endl;


    float i[] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    TRealVector<float, unsigned int>::vwrite(13,i);
    TRealVector<float, unsigned int>::vrev(13,i);
    TRealVector<float, unsigned int>::vwrite(13,i);
}

#endif
