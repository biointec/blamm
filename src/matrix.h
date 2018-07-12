/***************************************************************************
 *   Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)            *
 *   This file is part of Blamm                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cassert>
#include <vector>
#include <config.h>
#include <cstring>

#ifdef HAVE_MKL
        #include "mkl.h"
#endif

#ifdef HAVE_CUDA
        #include <cublas_v2.h>
#endif

// ============================================================================
// BLAS SINGLE/DOUBLE PRECISION FUNCTION PROTOTYPES
// ============================================================================

#define sgemm_f77 F77_FUNC (sgemm, SGEMM)

// general matrix-matrix multiplication
extern "C" void sgemm_f77(const char* transA, const char* transB,
                          const int* m, const int* n, const int* k,
                          const float* alpha,
                          const float* A, const int* LDA,
                          const float* B, const int* LDB,
                          const float* beta,
                          float *C, const int* LDC);

// ===========================================================================
// SGEMM BATCH PARAMETERS
// ===========================================================================

class SgemmBatchParams
{
public:
        int groupCount;
        int *m, *n, *k, *LDA, *LDB, *LDC;
        float *alpha, *beta;
        float **A_array, **B_array, **C_array;
#ifdef HAVE_MKL
        CBLAS_TRANSPOSE *trans;
        int *groupSize;
#endif

        /**
         * Constructor
         * @param groupCount Number of BLAS sgemm operations to perform
         */
        SgemmBatchParams(int groupCount) : groupCount(groupCount) {
                m = new int[groupCount];
                n = new int[groupCount];
                k = new int[groupCount];
                LDA = new int[groupCount];
                LDB = new int[groupCount];
                LDC = new int[groupCount];
                alpha = new float[groupCount];
                beta = new float[groupCount];
                A_array = new float*[groupCount];
                B_array = new float*[groupCount];
                C_array = new float*[groupCount];

#ifdef HAVE_MKL
                trans = new CBLAS_TRANSPOSE[groupCount];
                groupSize = new int[groupCount];
                for (int i = 0; i < groupCount; i++) {
                        trans[i] = CblasNoTrans;
                        groupSize[i] = 1;
                }
#endif
        }

        /**
         * Destructor
         */
        ~SgemmBatchParams() {
                delete [] m;
                delete [] n;
                delete [] k;
                delete [] LDA;
                delete [] LDB;
                delete [] LDC;
                delete [] alpha;
                delete [] beta;
                delete [] A_array;
                delete [] B_array;
                delete [] C_array;
#ifdef HAVE_MKL
                delete [] trans;
                delete [] groupSize;
#endif
        }
};

// ===========================================================================
// MATRIX CLASS (Column major storage)
// ===========================================================================

class Matrix
{
private:
        size_t rows;    // number of rows
        size_t cols;    // number of columns
        float *data;    // actual storage for the elements

        /**
         * Allocate memory for data
         */
        void allocateMemory() {
#ifdef HAVE_MKL
                data = (float*)mkl_malloc(rows*cols*sizeof(float), 64);
#else
                data = new float[rows*cols];
#endif
        }

        /**
         * Free memory for data
         */
        void freeMemory() {
                if (data == NULL)
                        return;
#ifdef HAVE_MKL
                mkl_free(data);
#else
                delete [] data;
#endif
        }

public:
        /**
         * Default constructor
         */
        Matrix() : rows(0), cols(0), data(NULL) {}

        /**
         * Create a nRows x nCols matrix
         * @param nRows Number of rows in the matrix
         * @param nCols Number of columns in the matrix
         */
        Matrix(size_t nRows, size_t nCols) : rows(nRows), cols(nCols) {
                assert(rows > 0);
                assert(cols > 0);
                allocateMemory();
        }

        /**
         * Create a nRows x nCols matrix and initialize it
         * @param nRows Number of rows in the matrix
         * @param nCols Number of columns in the matrix
         * @param el Initializer object
         */
        Matrix(size_t nRows, size_t nCols, float el) :
                Matrix(nRows, nCols) {
                fill(el);
        }

        /**
         * Forbid copy constructor
         * @param M Matrix to copy
         */
        Matrix(const Matrix& M) = delete;

        /**
         * Matrix destructor
         */
        ~Matrix() {
                freeMemory();
        }

        /**
         * Resize the matrix and initialize it
         * @param nRows Number of rows in the matrix
         * @param nCols Number of columns in the matrix
         * @param el Initializer object
         */
        void resize(size_t nRows, size_t nCols, float el) {
                assert(nRows > 0);
                assert(nCols > 0);

                // delete existing data
                freeMemory();

                // create memory for new matrix
                rows = nRows;
                cols = nCols;
                allocateMemory();
                fill(el);
        }

        /**
         * Set the matrix elements to zero fast
         */
        void setZero() {
                memset(data, 0, sizeof(float)*rows*cols);
        }

        /**
         * Fill the matrix with a certain element
         * @param el Element to fill the matrix with
         */
        void fill(const float& el) {
                for (size_t i = 0; i < rows*cols; i++)
                        data[i] = el;
        }

        /**
         * Retrieve the number of rows in the matrix
         * @return Number of rows
         */
        size_t nRows() const {
                return rows;
        }

        /**
         * Retrieve the number of columns in the matrix
         * @return Number of columns
         */
        size_t nCols() const {
                return cols;
        }

        /**
         * Overloaded parentheses to access/modify elements
         * @param row Row specification
         * @param col Column specification
         * @return Reference to element at specified position
         */
        float& operator()(size_t row, size_t col) {
                return data[col*rows+row];
        }

        /**
         * Overloaded parentheses to access/modify elements
         * @param row Row specification
         * @param col Column specification
         * @return Const-reference to element at specified position
         */
        const float& operator()(size_t row, size_t col) const {
                return data[col*rows+row];
        }

        /**
         * Get the data pointer
         * @return The data pointer
         */
        float* getData() const {
                return data;
        }

        /**
         * Print matrix elements to the output stream
         * @param os Output stream to add to
         * @param M Matrix to print
         * @return Output stream with the matrix elements
         */
        friend std::ostream& (operator<<)(std::ostream& os, const Matrix& M);

        /**
         * Print the sequence imposed by the matrix
         * @param overlap Number of overlapping nucleotides
         */
        void printSequence(size_t overlap) const;

        /**
         * Perform batch matrix-matrix multiplications
         * @param p sgemm batch paramters
         */
        static void sgemm_batch(const SgemmBatchParams& p)
        {
#ifdef HAVE_MKL
                cblas_sgemm_batch(CblasColMajor, p.trans, p.trans, p.m, p.n,
                                  p.k, p.alpha, (const float**)p.A_array, p.LDA,
                                  (const float**)p.B_array, p.LDB, p.beta,
                                  p.C_array, p.LDC, p.groupCount, p.groupSize);
#else
                for (int i = 0; i < p.groupCount; i++)
                        sgemm_f77("N", "N", &p.m[i], &p.n[i], &p.k[i],
                                  &p.alpha[i], p.A_array[i], &p.LDA[i],
                                  p.B_array[i], &p.LDB[i], &p.beta[i],
                                  p.C_array[i], &p.LDC[i]);
#endif
        }

#ifdef HAVE_CUDA
        /**
         * Perform batch matrix-matrix multiplications
         * @param handle cublas handle
         * @param p sgemm batch paramters
         */
        static void sgemm_batch_cuda(cublasHandle_t handle,
                                     const SgemmBatchParams& p)
        {
                for (int i = 0; i < p.groupCount; i++)
                        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, p.m[i],
                                    p.n[i], p.k[i], &p.alpha[i], p.A_array[i],
                                    p.LDA[i], p.B_array[i], p.LDB[i], &p.beta[i],
                                    p.C_array[i], p.LDC[i]);

        }
#endif
};

#endif
