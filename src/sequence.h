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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <deque>
#include <climits>
#include <map>
#include <set>
#include <fstream>
#include <mutex>

#include "matrix.h"
#include "motif.h"

// ============================================================================
// FUNCTION PROTOTYPE
// ============================================================================

class FastaBatch;

// ============================================================================
// SEQUENCE POSITION
// ============================================================================

class SeqPos {

private:
        size_t seqIdx;  // sequence index
        size_t seqPos;  // position in the sequence

public:
        /**
         * Default constructor
         */
        SeqPos() : seqIdx(0), seqPos(0) {}

        /**
         * Default constructor
         * @param seqIdx_ Sequence index
         * @param seqPos_ Sequence position
         */
        SeqPos(size_t seqIdx_, size_t seqPos_) :
                seqIdx(seqIdx_), seqPos(seqPos_) {}

        /**
         * Get the sequence index
         * @return The sequence index
         */
        size_t getSeqIndex() const {
                return seqIdx;
        }

        /**
         * Get the sequence position
         * @return The sequence position
         */
        size_t getSeqPos() const {
                return seqPos;
        }

        /**
         * Increment the sequence position by an offset
         * @param offset Offset
         */
        void incSeqPos(size_t offset = 1) {
                seqPos += offset;
        }
};

// ============================================================================
// SEQUENCE BLOCK
// ============================================================================

class SeqBlock {

private:
        /**
         * Check whether a SeqPos marker is contiguous with the final sequence
         * in the SeqBlock (= same contigous sequence)
         * @param sp SeqPos marker
         * @return true or false
         */
        bool isContiguous(const SeqPos& sp) const;

        std::string block;              // block of one or more sequences
        std::map<size_t, SeqPos> block2seq;     // block to sequence map

public:
        /**
         * Given a block position, get the corresponding sequence position
         * @param blockPos Block position
         * @return The sequence position
         */
        SeqPos getSeqPos(size_t blockPos) const;

        /**
         * Given a block position, get the remaining length within the sequence
         * @param blockPos Block position
         * @return The remaining length within this sequence
         */
        size_t getRemainingSeqLen(size_t blockPos) const;

        /**
         *
         *
         */
        std::string substr(size_t startPos, size_t len) {
                return block.substr(startPos, len);
        }

        /**
         * Append a new sequence to the block
         * @param data Sequence content
         * @param sp Sequence position
         */
        void append(const std::string& data, const SeqPos& sp);

        /**
         * Clear the sequence block
         */
        void clear() {
                block.clear();
                block2seq.clear();
        }

        /**
         * Is the sequence block empty()
         * @return True if empty
         */
        bool empty() const {
                return block.empty();
        }

        /**
         * Get the size of the sequence block
         * @return The size of the sequence block
         */
        size_t size() const {
                return block.size();
        }

        /**
         * Get a specific character
         * @param blockPos Block position
         */
        char operator[](size_t blockPos) {
                return block[blockPos];
        }

        /**
         * Copy the suffix of this block into a new block
         * @param newBlock New block (output)
         * @param startPos Start position of the suffix
         */
        void getSuffixBlock(SeqBlock& newBlock, size_t startPos);

        /**
         * operator<< overloading
         * @param os Output stream
         * @param sb Sequence block
         */
        friend std::ostream& operator<< (std::ostream& os, const SeqBlock& sb);
};

// ============================================================================
// FASTA BATCH HANDLER
// ============================================================================

class FastaBatch {

private:
        /**
         * Close the current input file stream and move to the next one
         * @return False if no more files are left, true otherwise
         */
        bool moveToNextFile();

        // following variables are modified as side effect of moveToNextFile()
        std::ifstream _ifs;                     // current input file stream
        size_t _currFileIdx;                    // current file index

        // following variables are modified as side effect of getNextLine()
        std::vector<std::string> _seqNames;     // sequence names
        size_t _currSeqLen;                     // current sequence length
        size_t _totSeqLen;                      // total sequence length

        /**
         * Get the next line of sequence content
         * @param line String where the line is stored
         * @param seqPos SeqPos structure
         * @return False if no more lines are left, true otherwise
         */
        bool getNextLine(std::string& line, SeqPos& seqPos);

        /**
         * Filter (and break) a line in valid ACTG sequences
         * @param line String with the input sequence content (unfiltered)
         * @param seqPos SeqPos structure
         * @param filtLine Filtered sequence lines (output)
         * @param filtSeqPos Filtered sequence positions (output)
         */
        void filterLine(const std::string& line, const SeqPos& seqPos,
                        std::deque<std::string>& filtLine,
                        std::deque<SeqPos>& filtSeqPos) const;

        /**
         * Append a sequence block with additional (raw) data
         * @param block Block of bulk sequence data
         * @param maxSize Maximum size of the block
         * @return True if at least one character was appended
         */
        bool appendNextBlock(SeqBlock& block, size_t maxSize);

        // following variables are modified as side effect of appendNextBlock()
        std::deque<std::string> _filtSeqv;
        std::deque<SeqPos> _filtSeqPosv;
        size_t _totFiltSeqLen;          // total filted sequence length
        size_t _maxFiltSeqLen;          // maximum filtered sequence length

        // following variables are modified as side effect of getNextOverlappingBlock()
        SeqBlock _nextBlock;
        std::mutex m;

        std::vector<std::string> seqFiles;      // sequence files in fasta batch

public:
        /**
         * Default constructor
         * @param seqFiles List of fasta files to read from
         * @param maxFiltSeqLen Maximum filtered sequence length
         */
        FastaBatch(const std::vector<std::string>& seqFiles,
                   size_t maxFiltSeqLen = std::numeric_limits<size_t>::max()) :
                   _currFileIdx(0), _currSeqLen(0), _totSeqLen(0),
                   _totFiltSeqLen(0), _maxFiltSeqLen(maxFiltSeqLen),
                   seqFiles(seqFiles) {}

        /**
         * Get the progress in percentage
         * @return Progress in percentage
         */
        float getProgressPerc() const {
                return (float)_totFiltSeqLen / (float)_maxFiltSeqLen;
        }

        /**
         * Get the sequence names
         */
        const std::vector<std::string>& getSeqNames() const {
                return _seqNames;
        }

        /**
         * Given a sequence index, get the sequence name
         * @param seqIdx Sequence index
         * @return The sequence name
         */
        std::string getSeqName(size_t seqIdx) const {
                return _seqNames[seqIdx];
        }

        /**
         * Get the number of sequence files in this batch
         * @return The number of sequence files in this batch
         */
        size_t size() const {
                return seqFiles.size();
        }

        /**
         * Get the total number of sequences in all files
         * @return The total number of sequences
         */
        size_t getNumSequences() const {
                return _seqNames.size();
        }

        /**
         * Get the total sequence length
         * @return The total sequence length
         */
        size_t getTotalSeqLength() const {
                return _totSeqLen;
        }

        /**
         * Get a filtered line of sequence data
         * @param line Line of sequence data
         * @param seqPos Sequence position
         * @return True if the line is non-empty
         */
        bool getNextFilteredLine(std::string& line, SeqPos& seqPos);

        /**
         * Get filtered sequence data
         * @param block Block of sequence data (max size = payload + overlap)
         * @param payload Desired size of the content
         * @param overlap Overlap between consecutive blocks
         * @return True if the block is non-empty
         */
        bool getNextOverlappingBlock(SeqBlock& block, size_t payload,
                                     size_t overlap);
};

// ============================================================================
// SEQUENCE MATRIX
// ============================================================================

class SeqMatrix {

public:
//private:
        // matrix dimensions = 4(K + overlap) x numCol
        size_t K;                       // number of sequence columns
        size_t overlap;                 // number of overlap columns
        size_t numRow;                  // number of rows
        size_t numOccRow;               // number of occupied rows

        Matrix S;                // actual matrix
        SeqBlock block;                 // sequence block

public:
        /**
         * Default constructor
         * @param numRow Number of sequence columns
         * @param K Number of sequence columns
         * @param overlap Number of overlap cols
         */
        SeqMatrix(size_t numRow, size_t K, size_t overlap) :
                K(K), overlap(overlap), numRow(numRow), numOccRow(0),
                S(numRow, 4*(K+overlap)) {}

        /**
         * Fill the next sequence matrix
         * @param bf Fasta batch reader
         * @return True if more data is present
         */
        bool getNextSeqMatrix(FastaBatch& bf);

        /**
         * Given a matrix position (row, col), get the remaining sequence length
         * @param row row index
         * @param col column index
         * @return The remaining length within this sequence
         */
        size_t getRemainingSeqLen(size_t row, size_t col) const {
                size_t blockPos = K*row + col;
                if (blockPos >= block.size())
                        return 0;
                return block.getRemainingSeqLen(blockPos);
        }

        float* getData() const {
                return S.getData();
        }

        /**
         * Given a matrix position (row, col), get the corresponding sequence position
         * @param row row index
         * @param col column index
         * @return The sequence position at (row, col)
         */
        SeqPos getSeqPos(size_t row, size_t col) const {
                return block.getSeqPos(K*row + col);
        }

        /**
         * Get the number of occupied columns in the matrix
         * @return The number of occupied columns
         */
        size_t getNumOccRow() const {
                return numOccRow;
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param sm Sequence matrix
         */
        friend std::ostream& operator<< (std::ostream& os, const SeqMatrix& sm);
};

#endif
