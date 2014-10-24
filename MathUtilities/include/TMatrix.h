/*
 * TMatrix is a templatized class that is to be used for typical matrix
 * operations.  The input to the template is the data type, usually a float or
 * double.
 *
 * Matrix algorithms taken from:
 *
 * 1.  Numerical Linear Algebra, Lloyd N. Trefethen and David Bau, SIAM 1997.
 * 2.  Introductory Linear Algebra with Applications, Bernard Kolman,
 *     Macmillan 1993.
 *
 * 2013-05-26   Geoff Hayes     Initial Release.
 */
#ifndef TMATRIX_H_
#define TMATRIX_H_

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <new>
#include <stdexcept>

//! Templatized class to be used for typical matrix operations.

/*!
  TMatrix is a templatized class that is to be used for typical matrix
  operations.  The input to the template is the data type, usually a float or
  double.

  Matrix algorithms taken from:

   1.  Numerical Linear Algebra, Lloyd N. Trefethen and David Bau, SIAM 1997.
   2.  Introductory Linear Algebra with Applications, Bernard Kolman,
       Macmillan 1993.
 */
template <class T>
class TMatrix
{
    public:

        /**
         * TMatrix class constructor.  Creates the matrix given the input
         * number of rows and columns.
         *
         * @param   rows  A positive number of rows for the matrix
         *                (default is zero).
         * @param   cols  A positive number of columns for the matrix
         *                (default is zero).
         *
         * @throw   std::bad_alloc if memory cannot be allocated to the
         *          matrix.
         */
        TMatrix(unsigned int rows = 0, unsigned int cols = 0) :
            _data(0), _numRows(rows), _numCols(cols)
        {
            // allocate memory for the matrix
            if (_numRows>0 && _numCols>0)
            {
                _data = new T[_numRows*_numCols];

                // throw an error if the memory couldn't be allocated
                if (!_data)
                {
                    _numRows = 0;
                    _numCols = 0;
                    throw std::bad_alloc();
                }
            }
            else
            {
                _numRows = 0;
                _numCols = 0;
            }

            // set the memory to zero
            memset(_data, 0, sizeof(T)*_numRows*_numCols);

        }

        /**
         * TMatrix class copy constructor.  Performs a deep copy of
         * the input matrix.
         *
         * @param   copy  The matrix to copy.
         *
         * @throw   std::bad_alloc if memory cannot be allocated to the
         *          matrix.
         *
         * @note    Invokes the TMatrix assignment operator.
         */
        TMatrix(const TMatrix& copy) :
            _data(0), _numRows(0), _numCols(0)
        {
            *this = copy;
        }

        /**
         * TMatrix class destructor.  Frees all memory allocated to
         * the matrix.
         */
        ~TMatrix()
        {
            if (_data)
            {
                delete [] _data;
                _data = 0;
            }
        }

        /**
         * Copies the contents of an array into the matrix given the input
         * number of rows and columns.
         *
         * @param   numRows   The number of rows in the data array.
         * @param   numCols   The number of columns in the data array.
         * @param   data      The array data to be copied into the matrix.
         *
         * @throw   std::bad_alloc if memory cannot be allocated to the matrix.
         */
        void copy (unsigned int numRows, unsigned int numCols, const T* data)
        {
            // allocated memory only if there is a difference in the array size
            // from the matrix size
            if (!_data || _numRows*_numCols != numRows*numCols)
            {
                delete [] _data;
                _data = 0;
                _data = new T(numRows*numCols);

                if (!_data)
                {
                    throw std::bad_alloc();
                }
            }

            memcpy(_data, data, sizeof(T)*numRows*numCols);

            _numRows = numRows;
            _numCols = numCols;
        }

        /**
         * Returns the integer number of rows of the matrix.
         */
        unsigned int getRows() const {return _numRows;}

        /**
         * Returns the integer number of columns of the matrix.
         */
        unsigned int getCols() const {return _numCols;}

        /**
         * Returns the element given by the row and column indices.
         *
         * @param   row   The positive matrix row index.
         * @param   col   The positive matrix column index (default is one).
         *
         * @throw   std::logic_error if the input row and/or column indices are
         *          invalid, or if no memory has been allocated to the matrix.
         */
        T& operator()(const unsigned int row, unsigned int col=1)
        {
            if (_data)
            {
                if (row>_numRows || row==0 || col>_numCols || col==0)
                {
                    throw std::logic_error(
                            "TMatrix::() - Invalid matrix accessor pair");
                }
                else
                {
                    return _data[(row-1)*_numCols + (col-1)];
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::() - Matrix has zero dimension");
            }
        }

        /**
         * TMatrix class assignment operator.  Performs a deep copy of the
         * input matrix.  Guards against self-assignment.
         *
         * @param   copy   The matrix to copy.
         *
         * @throw   std::bad_alloc if memory cannot be allocated to the
         *          matrix.
         */
        TMatrix& operator=(const TMatrix& copy)
        {
            if (this != &copy)
            {
                if (!_data ||
                        _numRows*_numCols != (copy._numRows*copy._numCols))
                {
                    // free any allocated memory
                    delete [] _data;
                    _data = 0;

                    // allocate memory given the new number of rows and columns
                    _data = new T[copy._numRows*copy._numCols];

                    if (!_data)
                    {
                        throw std::bad_alloc();
                    }
                }

                _numRows = copy._numRows;
                _numCols = copy._numCols;
                memcpy(_data, copy._data, sizeof(T)*_numRows*_numCols);
            }

            return *this;
        }

        /**
         * Resizes the matrix row and columns.  Deletes and re-allocates memory
         * if the new dimensions are not compatible with the existing ones.  All
         * matrix elements are cleared (set to zero).
         *
         * @param   rows  A positive number of rows for the matrix
         *                (default is zero).
         * @param   cols  A positive number of columns for the matrix
         *                (default is zero).
         *
         * @throw   std::bad_alloc if memory cannot be allocated to the
         *          matrix.
         */
        void resize(unsigned int rows=0, unsigned int cols=0)
        {
            if ((rows*cols) != (_numRows*_numCols))
            {
                // free any allocated memory
                delete [] _data;
                _data = 0;

                // allocate memory for the matrix
                if (rows>0 && cols>0)
                {
                    _data = new T[rows*cols];

                    // throw an error if the memory couldn't be allocated
                    if (!_data)
                    {
                        _numRows = 0;
                        _numCols = 0;
                        throw std::bad_alloc();
                    }
                }
                else
                {
                    rows = 0;
                    cols = 0;
                }
            }

            _numRows = rows;
            _numCols = cols;

            // set the memory/elements to zero
            memset(_data, 0, sizeof(T)*_numRows*_numCols);
        }

        /**
         * Writes the matrix contents to the console using std::cout.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix.
         */
        void write() const
        {
            if (_data)
            {
                std::cout << std::fixed;

                for (unsigned int i=0; i<_numRows; ++i)
                {
                    for (unsigned int j=0; j<_numCols; ++j)
                    {
                        std::cout << std::setprecision(8) << _data[i*_numCols + j] << "\t";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::write - Matrix has zero dimension");
            }
        }

        /**
         * Assigns a constant value to all elements in the matrix
         *
         * @param   value   Constant value to assign.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix.
         */
        TMatrix& operator=(const T value)
        {
            if (_data)
            {
                T* ptr = _data;
                for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                {
                    *ptr++ = value;
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Adds a constant value to all elements in the matrix
         *
         * @param   value  Constant value to add.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix.
         */
        TMatrix& operator+=(const T value)
        {
            if (_data)
            {
                T* ptr = _data;
                for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                {
                    *ptr++ += value;
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::+= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Subtracts a constant value from all elements in the matrix
         *
         * @param   value   Constant value to subtract.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix.
         */
        TMatrix& operator-=(const T value)
        {
            if (_data)
            {
                T* ptr = _data;
                for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                {
                    *ptr++ -= value;
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::-= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Multiplies a constant value against all elements in the matrix
         *
         * @param   value   Constant value to multiply.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix.
         */
        TMatrix& operator*=(const T value)
        {
            if (_data)
            {
                T* ptr = _data;
                for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                {
                    *ptr++ *= value;
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::*= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Divides a constant value against all elements in the matrix
         *
         * @param   value   Constant value to divide.
         *
         * @throw   std::logic_error if no memory has been allocated to the
         *          matrix, or if input value is zero.
         */
        TMatrix& operator/=(const T value)
        {
            if (_data)
            {
                if (value != static_cast<T>(0))
                {
                    T* ptr = _data;
                    for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                    {
                        *ptr++ /= value;
                    }
                }
                else
                {
                    throw std::logic_error(
                            "TMatrix::/= - Matrix division by zero");
                }

            }
            else
            {
                throw std::logic_error(
                        "TMatrix::/= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Performs element-by-element addition of two identically sized
         * matrices.
         *
         * @param   mtxToAdd   The matrix to add.
         *
         * @throw   std::logic_error if the two matrices are not identically
         *          sized (i.e. same number of rows and columns) or if no memory
         *          has been allocated to the "self" matrix.
         *
         * @return  Returns a reference to the updated (self) matrix.
         */
        TMatrix& operator+=(const TMatrix& mtxToAdd)
        {
            if (_data )
            {
                if (_numRows==mtxToAdd._numRows && _numCols==mtxToAdd._numCols)
                {
                    T* ptrA = _data;
                    T* ptrB = mtxToAdd._data;
                    for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                    {
                        *ptrA++ += *ptrB++;
                    }
                }
                else
                {
                    throw std::logic_error(
                            "TMatrix::+= - incompatible matrices");
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::+= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Performs element-by-element addition of two identically sized
         * matrices.
         *
         * @param   mtxToAdd   The matrix to add.
         *
         * @throw   std::logic_error if the two matrices are not identically
         *          sized (i.e. same number of rows and columns) or if no memory
         *          has been allocated to the "self" matrix.
         *
         *          std::bad_alloc if memory cannot be assigned to the matrix copy.
         *
         * @return  Returns a copy of the updated (self) matrix.
         */
        TMatrix operator+(const TMatrix& mtxToAdd)
        {
            if (_data )
            {
                TMatrix temp(*this);
                temp += mtxToAdd;
                return temp;
            }
            else
            {
                throw std::logic_error("TMatrix::+ - Matrix has zero dimension");
            }
        }

        /**
         * Performs element-by-element subtraction of two identically sized
         * matrices.
         *
         * @param   mtxToSub   The matrix to subtract.
         *
         * @throw   std::logic_error if the two matrices are not identically
         *          sized (i.e. same number of rows and columns) or if no memory
         *          has been allocated to the "self" matrix.
         *
         * @return  Returns a reference to the updated (self) matrix.
         */
        TMatrix& operator-=(const TMatrix& mtxToSub)
        {
            if (_data)
            {
                if (_numRows==mtxToSub._numRows && _numCols==mtxToSub._numCols)
                {
                    T* ptrA = _data;
                    T* ptrB = mtxToSub._data;
                    for (unsigned int i=0; i<(_numRows*_numCols); ++i)
                    {
                        *ptrA++ -= *ptrB++;
                    }
                }
                else
                {
                    throw std::logic_error(
                            "TMatrix::-= - incompatible matrices");
                }
            }
            else
            {
                throw std::logic_error("TMatrix::-= - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Performs element-by-element subtraction of two identically sized
         * matrices.
         *
         * @param   mtxToSub   The matrix to subtract.
         *
         * @throw   std::logic_error if the two matrices are not identically
         *          sized (i.e. same number of rows and columns) or if no memory
         *          has been allocated to the "self" matrix.
         *
         *          std::bad_alloc if memory cannot be assigned to the matrix copy.
         *
         * @return  Returns a copy of the updated (self) matrix.
         */
        TMatrix operator-(const TMatrix& mtxToSub)
        {
            if (_data )
            {
                TMatrix temp = *this;
                temp -= mtxToSub;
                return temp;
            }
            else
            {
                throw std::logic_error("TMatrix::- - Matrix has zero dimension");
            }
        }

        /**
         * Performs matrix multiplication of two matrices provided that they
         * are of compatible dimensions (i.e. the number of columns in the "left"
         * matrix is identical to the number of rows in the "right" matrix; left
         * and right referring to the position of the matrix relative to the
         * operator).
         *
         * @param   mtxToMult   The matrix to multiply.
         *
         * @throw   std::logic_error if the two matrices are not compatibly
         *          sized for matrix multiplication or if no memory
         *          has been allocated to the "self" matrix.
         *
         * @return  Returns a reference to the updated (self) matrix.
         */
        TMatrix& operator*=(const TMatrix& mtxToMult)
        {
            if (_data)
            {
                if (_numCols==mtxToMult._numRows)
                {
                    TMatrix temp(_numRows, mtxToMult._numCols);
                    T* ptrC = temp._data;

                    for (unsigned int i=0; i<temp._numRows; ++i)
                    {
                        for (unsigned int j=0; j<temp._numCols; ++j)
                        {
                            T* ptrA = &_data[i*_numCols];     // ith row
                            T* ptrB = &mtxToMult._data[j];    // jth col
                            T  sum = static_cast<T>(0);

                            for (unsigned int k=0; k<_numCols; ++k)
                            {
                                sum += (*ptrA++)*(*ptrB);
                                ptrB += mtxToMult._numCols;
                            }

                            *ptrC++ = sum;
                        }
                    }

                    // copy the temp matrix
                    *this = temp;
                }
                else
                {
                    throw std::logic_error(
                            "TMatrix::*= - incompatible matrices");
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::*+ - Matrix has zero dimension");
            }

            return *this;
        }

        /**
         * Performs matrix multiplication of two matrices provided that they
         * are of compatible dimensions (i.e. the number of columns in the "left"
         * matrix is identical to the number of rows in the "right" matrix; left
         * and right referring to the position of the matrix relative to the
         * operator).
         *
         * @param   mtxToMult   The matrix to multiply.
         *
         * @throw   std::logic_error if the two matrices are not compatibly
         *          sized for matrix multiplication or if no memory
         *          has been allocated to the "self" matrix.
         *
         * @return  Returns a copy of the updated (self) matrix.
         */
        TMatrix operator*(const TMatrix& mtxToMult)
        {
            if (_data )
            {
                TMatrix temp = *this;
                temp *= mtxToMult;
                return temp;
            }
            else
            {
                throw std::logic_error("TMatrix::* - Matrix has zero dimension");
            }
        }

        /**
         * Calculates the transpose of the matrix.
         *
         * @throw    std::bad_alloc if memory cannot be allocated to the
         *           transposed matrix.
         * @throw    std::logic_error if the matrix to transpose is of zero
         *           dimension (no memory allocated to it).
         * @return   Retuns a copy of the transposed matrix.
         */
        TMatrix t() const
        {
            if (_data )
            {
                TMatrix temp(_numCols, _numRows);

                for (unsigned int i=0; i<_numRows; ++i)
                {
                    for (unsigned int j=0; j<_numCols; ++j)
                    {
                        temp._data[j*_numRows+i] = _data[i*_numCols+j];
                    }
                }

                return temp;
            }
            else
            {
                throw std::logic_error("TMatrix::t - Matrix has zero dimension");
            }
        }

        /**
         * Returns a boolean indicating if the matrix is symmetric or not.  A
         * matrix is considered symmetric if it is square and row i is identical
         * to column i.
         *
         * @retval   true if the matrix is symmetric.
         * @retval   false if the matrix is not symmetric.
         */
        bool isSymmetric() const
        {
            bool status = (_numRows==_numCols) && _numRows>0;

            if (status)
            {
                for (unsigned int i=0; i<_numRows; ++i)
                {
                    for (unsigned int j=0; j<_numCols; ++j)
                    {
                        if (_data[i*_numCols+j] !=
                                _data[j*_numRows+i])
                        {
                            status = false;
                            break;
                        }
                    }
		    
                    if (!status)
                    {
                        break;
                    }
                }
            }

            return status;
        }

        /**
         * Sets the matrix as the identity matrix i.e. ones along the diagonal
         * and zeros elsewhere.
         *
         * @throw   std::logic_error if the matrix is not square.
         */
        void setAsIdentity()
        {
            if ((_numRows==_numCols) && _numRows>0)
            {
                memset(_data, 0, sizeof(T)*_numRows*_numCols);

                for (unsigned int i=0; i<_numRows; ++i)
                {
                    _data[i*_numRows+i] = static_cast<T>(1);
                }
            }
            else
            {
                throw std::logic_error(
                        "TMatrix::setAsIdentity - Matrix is not square.");
            }
        }

        /**
         * Caclulates the LU factorization of the given matrix where L is a unit
         * lower triangular matrix and U is an upper triangular matrix.  This
         * algorithm makes use of partial pivoting (is more stable) than that
         * without and so provides a permutation matrix P. i.e. PA = LU where A
         * is the self matrix.
         *
         * @param   L   Unit lower triangular matrix (i.e. ones along its
         *              diagonal).
         * @param   U   Upper triangular matrix.
         * @param   P   Permutation matrix indicating the permutated rows of the
         *              original matrix.
         * @param   numPerms   Number of permutations that have been invoked
         *                     against P.
         *
         * @retval  true if the LU factorization was successful.
         * @retval  false if the LU factorization was unsuccessful.
         *
         * @throw   std::bad_alloc if memory cannot be allocated to either L, U
         *          or P.
         *
         * @note    Assumes that the "self" matrix is square.
         */
        bool getLUPFactorization(TMatrix& L, TMatrix& U, TMatrix& P,
                                unsigned int* numPerms = 0) const
        {
            // assume valid for square matrices only
            bool status = (_numRows==_numCols) && _numRows>0;

            if (numPerms)
            {
                *numPerms = 0;
            }

            if (status)
            {
                U = *this;
                L.resize(_numRows,_numCols);
                L.setAsIdentity();
                P.resize(_numRows,_numCols);
                P.setAsIdentity();

                for (unsigned int k=0; k<_numRows-1; ++k)
                {
                    // find that row of U such that |u{ik}| is maximized
                    T            maxUik       = std::abs(U._data[k*_numRows + k]);
                    unsigned int maxAt        = k;
                    bool         doInterchange = false;

                    for (unsigned int i=k+1; i<_numRows; ++i)
                    {
                        if (std::abs(U._data[i*_numRows + k]) > maxUik)
                        {
                            maxUik = std::abs(U._data[i*_numRows + k]);
                            maxAt  = i;
                            doInterchange = true;
                        }
                    }

                    // interchange the two rows if necessary
                    if (doInterchange)
                    {
                        for (unsigned int i=k; i<_numRows; ++i)
                        {
                            T temp = U._data[k*_numRows+i];
                            U._data[k*_numRows+i] = U._data[maxAt*_numRows+i];
                            U._data[maxAt*_numRows+i] = temp;
                        }

                        if (k>0)
                        {
                            for (unsigned int i=0; i<=k-1; ++i)
                            {
                                T temp = L._data[k*_numRows+i];
                                L._data[k*_numRows+i] = L._data[maxAt*_numRows+i];
                                L._data[maxAt*_numRows+i] = temp;
                            }
                        }

                        for (unsigned int i=0; i<_numRows; ++i)
                        {
                            T temp = P._data[k*_numRows+i];
                            P._data[k*_numRows+i] = P._data[maxAt*_numRows+i];
                            P._data[maxAt*_numRows+i] = temp;
                        }

                        if (numPerms)
                        {
                            *numPerms = *numPerms+1;
                        }
                    }


                    for (unsigned int j=k+1; j<_numRows; ++j)
                    {
                        T denom = U._data[k*_numRows+k];

                        if (std::abs(denom)==0.0)
                        {
                            // no solution
                            status = false;
                            break;
                        }

                        T temp = U._data[j*_numRows+k] / denom;

                        L._data[j*_numRows+k] = temp;

                        for (unsigned int i=k;i<_numRows;++i)
                        {
                            U._data[j*_numRows+i] -= temp*U._data[k*_numRows+i];
                        }
                    }
                }
            }

            return status;
        }

        /**
         * Function returns the determinant of a matrix.
         *
         * @return   The determinant of the matrix; 0 may mean that no
         *           determinant exists.
         *
         * @throw    std::logic_error if the matrix is not square.
         */
        T det() const
        {
            T determ = 0;

            if (_numRows==_numCols)
            {
                // calculate the determinant using the LUP factorization
                unsigned int numPerms = 0;
                TMatrix L;
                TMatrix U;
                TMatrix P;

                const bool isValid = getLUPFactorization(L,U,P,&numPerms);

                if (isValid)
                {
                    // we have calculated P*A = L*U => A = P^-1*U*P
                    // => determ(A) = determ(P^-1*U*P)
                    // => determ(A) = determ(P^-1)*determ(L)*determ(U)
                    // as P is a permuation matrix, then determ(P^-1) = (-1)^n
                    // where n is the number of rows that have been permuted in
                    // P
                    // as L and U are triangular, the determinant of each is
                    // just the product of the elements along the diagonal
                    //
                    // but since L is unit lower triangular, then its
                    // determinant is just one and so can be ignored
                    //
                    // but since the determinant of the permutation matrix P is
                    // just (-1)^numPerms, then if numPerms is positive then its
                    // determinant is +1 (and so can be ignored) else it is
                    // -1 and so has to be applied.
                    T determU = U._data[0];

                    for (unsigned int i=1;i<_numRows;++i)
                    {
                        determU *= U._data[i*_numRows + i];
                    }

                    determ = determU;

                    // if the number of permuations is odd, then negate the
                    // calculated determinant
                    if (numPerms % 2 != 0)
                    {
                        determ = -determ;
                    }
                }
            }
            else
            {
                throw std::logic_error("TMatrix::det - Matrix is not square.");
            }

            return determ;
        }

        /**
         * Calculates the inverse of the assumed square matrix using the LUP
         * factorization.  Note that it is assumed that the diagonal elements
         * along the L and U matrices are non-zero.  (If any should happen to
         * be zero, then the matrix is considered invertible.)
         *
         * @throw std::logic_error if the matrix is not square, no invertible
         *                         matrix exists (i.e. matrix is singular) or
         *                         no LUP factorization exists.
         *        std::bad_alloc   if memory for the matrix inverse cannot be
         *                         allocated (or for any temporary matrix)
         * @return                 A copy of the inverse matrix.
         */
        TMatrix i() const
        {
            if (_numRows==_numCols)
            {
                // calculate the LUP factorization
                TMatrix L;
                TMatrix U;
                TMatrix P;

                const bool status = getLUPFactorization(L,U,P);

                if (status)
                {
                    // create the inverse matrix
                    TMatrix inverse(_numRows, _numCols);

                    // create a vector to solve against
                    TMatrix b(_numRows, 1);

                    // temp vectors
                    TMatrix y(_numRows, 1);
                    TMatrix x(_numRows, 1);

                    // to solve the matrix inverse problem we will use forward
                    // and backwards substitution to solve for it:
                    // inv(A) = inv(L)*inv(U)*P
                    // b will be the vector that will be used for each column in
                    // the identity matrix
                    for (unsigned int i=0; i<_numRows; ++i)
                    {
                        // set b as the ith column of the identity matrix
                        memset(b._data, 0, sizeof(T)*_numRows);
                        b._data[i] = static_cast<T>(1);

                        // compute P*b
                        const TMatrix Pb = P*b;

                        // use forward substitution to solve for Ly = Pb
                        for (unsigned int j=0; j<_numRows; ++j)
                        {
                            const T denom = L._data[j*_numRows+j];

                            if (denom==static_cast<T>(0))
                            {
                                throw std::logic_error(
                                        "TMatrix::i - matrix is not invertible");
                            }

                            y._data[j] = Pb._data[j]/denom;

                            for (unsigned int k=0; k<j; ++k)
                            {
                                y._data[j] -= L._data[j*_numRows+k]*y._data[k]/denom;
                            }
                        }

                        // use backward substitution to solve for Ux = y
                        for (int j=_numRows-1; j>=0; --j)
                        {
                            const T denom = U._data[j*_numRows+j];

                            if (denom==static_cast<T>(0))
                            {
                                throw std::logic_error(
                                        "TMatrix::i - matrix is not invertible");
                            }

                            x._data[j] = y._data[j]/denom;

                            for (unsigned int k=j+1; k<_numRows; ++k)
                            {
                                x._data[j] -=
                                        U._data[j*_numRows+k]*x._data[k]/denom;
                            }
                        }

                        // update the ith column of the inverse matrix with x
                        for (unsigned int j=0; j<_numRows;++j)
                        {
                            inverse._data[j*_numRows+i] = x._data[j];
                        }
                    }

                    return inverse;
                }
                else
                {
                    throw std::logic_error(
                            "TMatrix::i - LUP Factorization failed.");
                }
            }
            else
            {
                throw std::logic_error("TMatrix::i - Matrix is not square");
            }
        }

        /**
         * Caclulates the QR factorization of the given matrix using the
         * modified Gram-Schmidt algorithm: Q is an orthongonal matrix and R
         * is an upper triangular matrix.
         *
         * @param   Q   Orthongonal matrix.
         * @param   R   Upper triangular matrix.
         *
         * @retval  true if the QR factorization was successful.
         * @retval  false if the QR factorization was unsuccessful.
         *
         * @throw   std::bad_alloc if memory cannot be allocated to either Q or
         *          R.
         *
         * @note    Assumes that the "self" matrix is square.
         */
        bool getQRFactorization(TMatrix& Q, TMatrix& R) const
        {
            bool status = (_numRows==_numCols) && _numRows>0;

            if (status)
            {
                // still need to implement
            }


            return false;
        }

    private:

        //! Pointer to memory allocated for the matrix data "store".
        T* _data;

        //! Integer number of rows in the matrix.
        unsigned int _numRows;

        //! Integer number of columns in the matrix.
        unsigned int _numCols;
};




#endif /* TMATRIX_H_ */
