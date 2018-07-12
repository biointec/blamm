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

#include "matrix.h"

using namespace std;

// ===========================================================================
// MATRIX CLASS
// ===========================================================================

void Matrix::printSequence(size_t overlap) const
{
        cout << "Print matrix: " << nRows() << " x " << nCols() << endl;
        int K = nCols() / 4 - overlap;
        cout << "K: " << K << ", overlap: " << overlap << endl;

        for (size_t r = 0; r < nRows(); r++) {
                for (size_t c = 0; c < nCols(); c += 4) {
                        if ((*this)(r, c+0) == 1)
                                cout << "A\t";
                        if ((*this)(r, c+1) == 1)
                                cout << "C\t";
                        if ((*this)(r, c+2) == 1)
                                cout << "G\t";
                        if ((*this)(r, c+3) == 1)
                                cout << "T\t";
                }
        }

        cout << endl;
}

/**
 * Specialized print matrix elements to the output stream for floats
 * @param os Output stream to add to
 * @param M Matrix to print
 * @return Output stream with the matrix elements
 */
std::ostream& operator<<(std::ostream& os, const Matrix& M)
{
        std::streamsize oldPrec = os.precision();
        os.precision(2);
        for (size_t r = 0; r < M.nRows(); r++) {
                for (size_t c = 0; c < (M.nCols()-1); c++)
                        os << M.data[c*M.nRows()+r] << "\t";
                os << M.data[(M.nCols()-1)*M.nRows()+r] << endl;
        }
        cout.precision(oldPrec);
        return os;
}
