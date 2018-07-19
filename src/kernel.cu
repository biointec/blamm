/***************************************************************************
 *   Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)            *
 *   This file is part of Blamm                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
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

__global__
void filterScore(float* R, int m, int n, float* threshold,
                 int* occIdx, float* occScore, int* nOcc)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x;  // row
        int j = blockIdx.y * blockDim.y + threadIdx.y;  // col
        int idx = j*m + i;      // column-major storage

        if ((i < m) && (j < n) && (R[idx] >= threshold[j])) {
                int resPos = atomicAdd(nOcc, 1);
                occScore[resPos] = R[idx];
                occIdx[resPos] = idx;
        }
}

void kernel_wrapper(float *d_R, int m, int n, float *d_threshold,
                    int* d_occIdx, float* d_occScore, int *d_nOcc)
{
        dim3 threadsPerBlock(32, 32);
        dim3 numBlocks((m + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (n + threadsPerBlock.y - 1) / threadsPerBlock.y);

        filterScore<<<numBlocks, threadsPerBlock>>>(d_R, m, n, d_threshold,
                                                    d_occIdx, d_occScore, d_nOcc);
}
