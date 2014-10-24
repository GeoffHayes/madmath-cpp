/*
 * Class to represent a complex number given real and imaginary parts.
 *
 * 2013-05-31   Geoff Hayes     Initial Release.
 */

#ifndef TCOMPLEXNUMBER_H_
#define TCOMPLEXNUMBER_H_

#include <cmath>
#include <iomanip>
#include <iostream>

//! Templatized class to be used for float or double data type complex numbers.

/*!
  TComplexNumber is a templatized class that is to be used to represent complex
  numbers.
 */
template <class T>
class TComplexNumber
{
    public:

        /**
         * TComplexNumber class constructor.
         *
         * @param   real   The real part of the complex number (default to
         *                 zero).
         * @param   imag   The imaginary part of the complex number
         *                 (default to zero).
         */
        TComplexNumber(T real = static_cast<T>(0), T imag = static_cast<T>(0))
        : _real(real),
          _imag(imag)
        {
            // intentionally left blank
        }

        /**
         * TComplexNumber class copy constructor.
         *
         * @param   cmplxNum  The complex number to copy.
         */
        TComplexNumber(const TComplexNumber& cmplxNum)
        : _real(cmplxNum._real),
          _imag(cmplxNum._imag)
        {
            // intentionally left blank
        }

        /**
         * TComplexNumber class destructor.
         */
        ~TComplexNumber()
        {
            // intentionally left blank
        }

        /**
         * Returns the real part of the complex number.
         *
         * @retval  Real part of the complex number.
         */
        T getReal() const {return _real;}

        /**
         * Returns the imaginary part of the complex number.
         *
         * @retval   Imaginary part of the complex number.
         */
        T getComplex() const {return _imag;}

        /**
         * Sets the real part of the complex number.
         *
         * @param   real   The new real part of the complex number.
         */
        void setReal(T real){_real = real;}

        /**
         * Sets the imaginary part of the complex number.
         *
         * @param   imag   The new imaginary part of the complex number
         */
        void setImag(T imag){_imag = imag;}

        /**
         * Complex number addition.
         *
         * @param   cmplxNum   The complex number to add.
         *
         * @retval  A copy of the updated complex number.
         */
        TComplexNumber operator+(const TComplexNumber& cmplxNum)
        {
            TComplexNumber temp(*this);
            temp += cmplxNum;

            return temp;
        }

        /**
         * Complex number addition assignment.
         *
         * @param   cmplxNum   The complex number to add.
         *
         * @retval  A reference to the updated complex number.
         */
        TComplexNumber& operator+=(const TComplexNumber& cmplxNum)
        {
            _real += cmplxNum._real;
            _imag += cmplxNum._imag;

            return (*this);
        }

        /**
         * Complex number subtraction.
         *
         * @param   cmplxNum   The complex number to subtract.
         *
         * @retval  A copy of the updated complex number.
         */
        TComplexNumber operator-(const TComplexNumber& cmplxNum)
        {
            TComplexNumber temp(*this);
            temp -= cmplxNum;

            return temp;
        }

        /**
         * Complex number subtraction assignment.
         *
         * @param   cmplxNum   The complex number to subtraction.
         *
         * @retval  A reference to the updated complex number.
         */
        TComplexNumber& operator-=(const TComplexNumber& cmplxNum)
        {
            _real -= cmplxNum._real;
            _imag -= cmplxNum._imag;

            return (*this);
        }

        /**
         * Complex number multiplication.
         *
         * @param   cmplxNum   The complex number to multiply by.
         *
         * @retval  A copy of the updated complex number.
         */
        TComplexNumber operator*(const TComplexNumber& cmplxNum)
        {
            TComplexNumber temp(*this);
            temp *= cmplxNum;

            return temp;
        }

        /**
         * Complex number multiplcation assignment.
         *
         * @param   cmplxNum   The complex number to multiply by.
         *
         * @retval  A reference to the updated complex number.
         */
        TComplexNumber& operator*=(const TComplexNumber& cmplxNum)
        {
            T prevReal = _real;

            _real = _real*cmplxNum._real - _imag*cmplxNum._imag;
            _imag = _imag*cmplxNum._real + prevReal*cmplxNum._imag;

            return (*this);
        }

        /**
         * Complex number division.
         *
         * @param   cmplxNum   The complex number to divide by.
         *
         * @retval  A copy of the updated complex number.
         */
        TComplexNumber operator/(const TComplexNumber& cmplxNum)
        {
            TComplexNumber temp(*this);
            temp /= cmplxNum;

            return temp;
        }

        /**
         * Complex number division assignment.
         *
         * @param   cmplxNum   The complex number to divide by.
         *
         * @retval  A reference to the updated complex number.
         */
        TComplexNumber& operator/=(const TComplexNumber& cmplxNum)
        {
            const T mag     = cmplxNum.abs();
            const T magSqrd = mag*mag;

            (*this) *= cmplxNum.conjugate();
            _real /= magSqrd;
            _imag /= magSqrd;

            return (*this);
        }

        /**
         * Returns the conjugate of the complex number.
         *
         * @retval   A copy to the complex number conjugate.
         */
        TComplexNumber conjugate() const
        {
            TComplexNumber temp(*this);
            temp._imag = -temp._imag;

            return temp;
        }

        /**
         * Returns the absolute value (or modulus or magnitude) of the complex
         * number.
         *
         * @retval   The absolute value of the complex number.
         */
        T abs() const
        {
            return sqrt(_real*_real + _imag*_imag);
        }

        /**
         * Writes the complex number to the console.
         */
        void write() const
        {
            std::cout << std::fixed;
            std::cout << std::setprecision(8)
                      << _real
                      << (_imag < static_cast<T>(0) ? "" : "+")
                      << _imag << "i";
            std::cout << std::endl;
        }

    private:

        //! The real part representation of the complex number.
        T _real;

        //! The imaginary part representation of the complex number.
        T _imag;
};

#endif /* TCOMPLEXNUMBER_H_ */
