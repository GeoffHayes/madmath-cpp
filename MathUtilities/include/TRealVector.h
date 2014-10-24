/*
 * Real vector math routines.
 *
 * 2013-06-01   Geoff Hayes     Initial Release.
 */

#ifndef TREALVECTOR_H_
#define TREALVECTOR_H_

#include <cmath>
#include <fenv.h>
#include <iomanip>
#include <iostream>

//! Templatized real (float or double) vector math routines where T is the data
//! type being used for each element of each vector, and U is used to represent
//! the number of elements in each vector (unsigned short, long, etc.).

/*!
  TRealVector is a templatized class that is to be used to provide static member
  functions to perform real vector routines.
 */
template <class T, class U> class TRealVector
{
    public:

        /**
         * Adds two sets of two real vectors and multiplies the two resulting
         * real vectors together.
         *
         * @param   y   Vector y, where y=u*v + w*x.
         * @param   u   Vector u, summed with v.
         * @param   v   Vector v, summed with u.
         * @param   w   Vector w, summed with x.
         * @param   x   Vector x, summed with w.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vaam(U n, T* y, const T* u, const T* v, const T* w, const T* x)
        {
            if (y && u && v && w && x && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ = (*u++ + *v++)*(*w++ + *x++);
                }
            }
        }

        /**
         * Computes the absolute value for each element in the vector.
         *
         * @param   y   Vector y, where y=abs(u).
         * @param   u   Vector u.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vabs(U n, T* y, const T* u)
        {
            if (y && u && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ = std::abs(*u++);
                }
            }
        }

        /**
         * Computes the arccosine for each element in the vector, in radians
         * within the interval of [0, pi].
         *
         * @param   y   Vector y, where y=acos(u).
         * @param   u   Vector u.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         * @warning Assumes that each element within the input vector u falls in
         *          the interval [-1,+1].  Element is cleared (set to zero) if
         *          this assumption is not satisfied.
         */
        static void vacos(U n, T* y, const T* u)
        {
            if (y && u && n)
            {
                for (U i=0; i<n; ++i)
                {
                    if (*u > static_cast<T>(1) || *u < static_cast<T>(-1))
                    {
                        // u is invalid so increment pointers and clear y
                        *y++ = static_cast<T>(0);
                        u++;
                    }
                    else
                    {
                        *y++ = std::acos(*u++);
                    }
                }
            }
        }

        /**
         * Element-by-element real vector addition as y[i] = u[i] + v[i].
         *
         * @param   y   Vector y to be the sum of u and v.
         * @param   u   Vector u.
         * @param   v   Vector v.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vadd(U n, T* y, const T* u, const T* v)
        {
            if (y && u && v && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ += *u++ + *v++;
                }
            }
        }

        /**
         * Truncates each element within the vector to its integer part and
         * outputs the result as a real.
         *
         * @param   y   Vector y, where y=floor(u).
         * @param   u   Vector u.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vaint(U n, T* y, const T* u)
        {
 /*           if (y && u && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++= std::trunc(*u++);
                  if (*u < static_cast<T>(0) && temp < *u)
                    {
                        // integer portion is one less than what it should be so
                        // need to add one (i.e. floor(-1.3) = -2.0.
                        temp += static_cast<T>(1);
                    }

                    *y++ = temp;
                    u++;
                }
            }*/
        }

        /**
         * Adds two real vectors together and multiplies the sum by a third.
         *
         * @param   y   Vector y, where y=(u+v)*w.
         * @param   u   Vector u, summed with v.
         * @param   v   Vector v, summed with u.
         * @param   w   Vector w, summed with x.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vam(U n, T* y, const T* u, const T* v, const T* w)
        {
            if (y && u && v && w && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ = (*u++ + *v++)*(*w++);
                }
            }
        }

        /**
         * Rounds the elements of the vector to their integer part and outputs
         * them as reals.  Elements are rounded down to the nearest integer.
         *
         * @param   y   Vector y, where y=floor(u).
         * @param   u   Vector u.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vanint(U n, T* y, const T* u)
        {
  /*          std::fesetround(FE_DOWNWARD);
            if (y && u && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ = std::nearbyint(*u++);
                }
            }*/
        }

        /**
         * Element-by-element real vector assignment as a[i] = alpha.
         *
         * @param   y      Vector y.
         * @param   alpha  Scalar value to be assigned to all elements of a.
         * @param   n      Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vset(U n, T* y, T alpha)
        {
            if (y && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ = alpha;
                }
            }
        }

        /**
         * Writes each element of vector to the console.
         *
         * @param   y      Vector y.
         * @param   n      Number of elements in each vector.
         *
         */
        static void vwrite(U n, const T* y)
        {
            if (y && n)
            {
                std::cout << std::fixed;

                for (U i=0; i<n; ++i)
                {
                    std::cout << "y[" << i << "]="
                              << std::setprecision(8) << *y++;
                    std::cout << std::endl;
                }
            }
            else
            {
                std::cout << "Vector y is empty" << std::endl;
            }
        }

        /**
         * Element-by-element real vector addition as y[i] = y[i] + u[i].
         *
         * @param   y   Vector y.
         * @param   u   Vector u.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector y is updated via this function.
         */
        static void vadd(U n, T* y, const T* u)
        {
            if (y && u && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *y++ += *u++;
                }
            }
        }

        /**
         * Element-by-element real vector subtraction as a[i] = a[i] - b[i].
         *
         * @param   a   Vector a.
         * @param   b   Vector b.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vsub(U n, T* a, const T* b)
        {
            if (a && b && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *a++ -= *b++;
                }
            }
        }

        /**
         * Element-by-element real vector multiplication as a[i] = a[i] * b[i].
         *
         * @param   a   Vector a.
         * @param   b   Vector b.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vmult(U n, T* a, const T* b)
        {
            if (a && b && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *a++ *= *b++;
                }
            }
        }

        /**
         * Element-by-element real vector division as a[i] = a[i] / b[i].
         *
         * @param   a   Vector a.
         * @param   b   Vector b.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         * @warning There is no check for a division by zero error.
         */
        static void vdiv(U n, T* a, const T* b)
        {
            if (a && b && n)
            {
                for (U i=0; i<n; ++i)
                {
                    *a++ /= *b++;
                }
            }
        }

        /**
         * Element-by-element real vector dot product as sum(a[i]*b[i]).
         *
         * @param   a   Vector a.
         * @param   b   Vector b.
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static T vdot(U n, T* a, const T* b)
        {
            T result = static_cast<T>(0);

            if (a && b && n)
            {
                for (U i=0; i<n; ++i)
                {
                    result += (*a++) * (*b++);
                }
            }

            return result;
        }

        /**
         * Element-by-element real vector linear combination of a single vector
         * as: a[i] = alpha*b[i]
         *
         * @param   a       Vector a.
         * @param   b       Vector b.
         * @param   alpha   Scalar applied to each element of vector b.
         * @param   n       Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vlincomb(U n, T* a, const T* b, T alpha)
        {
            if (a && b && n)
            {
                for (U i=0; i<n; ++i)
                {
                    (*a++) = alpha * (*b++);
                }
            }
        }

        /**
         * Element-by-element real vector linear combination of a two vectors
         * as: a[i] = alpha*b[i] + beta*c[i]
         *
         * @param   a       Vector a.
         * @param   b       Vector b.
         * @param   c       Vector c.
         * @param   alpha   Scalar applied to each element of vector b.
         * @param   beta    Scalar applied to each element of vector c.
         * @param   n       Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vlincomb(U n, T* a, const T* b,
                             const T* c, T alpha, T beta)
        {
            if (a && b && c && n)
            {
                for (U i=0; i<n; ++i)
                {
                    (*a++) = alpha * (*b++) + beta * (*c++);
                }
            }
        }

        /**
         * Element-by-element real vector linear combination of a three vectors
         * as: a[i] = alpha*b[i] + beta*c[i] + gamma*d[i]
         *
         * @param   a       Vector a.
         * @param   b       Vector b.
         * @param   c       Vector c.
         * @param   d       Vector d.
         * @param   alpha   Scalar applied to each element of vector b.
         * @param   beta    Scalar applied to each element of vector c.
         * @param   gamma   Scalar applied to each element of vector d.
         * @param   n       Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vlincomb(U n, T* a, const T* b,
                             const T* c, const T* d, T alpha, T beta, T gamma)
        {
            if (a && b && c && n)
            {
                for (U i=0; i<n; ++i)
                {
                    (*a++) = alpha * (*b++) + beta * (*c++) + gamma * (*d++);
                }
            }
        }

        /**
         * Determines the sum of all elements within the vector.
         *
         * @param   a   Vector a;
         * @param   n   Number of elements in each vector.
         */
        static T vsum(U n, const T* a)
        {
            T sum = static_cast<T>(0);

            if (a && n)
            {
                for (U i=0; i<n; ++i)
                {
                    sum += *a++;
                }
            }

            return sum;
        }

        /**
         * Determines the sum of squares of all elements within the vector.
         *
         * @param   a   Vector a;
         * @param   n   Number of elements in each vector.
         */
        static T vsumsq(U n, const T* a)
        {
            T sumsq = static_cast<T>(0);

            if (a && n)
            {
                for (U i=0; i<n; ++i)
                {
                    T temp = *a++;
                    sumsq += temp*temp;
                }
            }

            return sumsq;
        }

        /**
         * Determines the mean of the elements within the vector.
         *
         * @param   a   Vector a;
         * @param   n   Number of elements in each vector.
         */
        static T vmean(U n, const T* a)
        {
            T mean = static_cast<T>(0);

            if (a && n)
            {
                const T sum = TRealVector::vsum(n,a);
                mean = sum/n;
            }

            return mean;
        }

        /**
         * Determines the standard deviation of the elements within the vector.
         *
         * @param   a   Vector a;
         * @param   n   Number of elements in each vector.
         */
        static T vstddev(U n, const T* a)
        {
            T stddev = static_cast<T>(0);

            if (a && n>1)
            {
                T sum   = static_cast<T>(0);
                T sumsq = static_cast<T>(0);

                for (U i=0; i<n; ++i)
                {
                    T temp = *a++;
                    sum   += temp;
                    sumsq += temp*temp;
                }

                stddev = sqrt((sumsq - sum*sum/n))/(n-1);
            }

            return stddev;
        }

        /**
         * Reverses all the elements in the vector i.e. the a[0] = a[n], a[1] =
         * a[n-1], etc.
         *
         * @param   a   Vector a;
         * @param   n   Number of elements in each vector.
         *
         * @note    Vector a is updated via this function.
         */
        static void vrev(U n, T* a)
        {
            if (a && n)
            {
                T* b = a+(n-1);

                for (U i=0; i<n/2; ++i)
                {
                    T temp = *a;
                    *a++   = *b;
                    *b--   = temp;
                }
            }
        }
};

#endif /* TREALVECTOR_H_ */
