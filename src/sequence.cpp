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

#include <stdexcept>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#include "sequence.h"

using namespace std;

// ============================================================================
// SEQUENCE BLOCK
// ============================================================================

bool SeqBlock::isContiguous(const SeqPos& sp) const
{
        // if there is nothing to compare with: no
        if (block2seq.empty())
                return false;

        // if it's a difference sequence: no
        const SeqPos& last = block2seq.rbegin()->second;
        if (last.getSeqIndex() != sp.getSeqIndex())
                return false;

        // if it's continguous
        size_t blockDelta = block.size() - block2seq.rbegin()->first;
        size_t seqDelta = sp.getSeqPos() - last.getSeqPos();
        return (blockDelta == seqDelta);
}

SeqPos SeqBlock::getSeqPos(size_t blockPos) const
{
        // check the validity of the input parameter
        // assert (blockPos < block.size());

        // find the last marker less than or equal to blockPos
        auto it = block2seq.upper_bound(blockPos);
        assert(it != block2seq.begin());
        it--;   // now contains last marker <= blockPos

        // now return the sequence position
        const SeqPos& last = it->second;
        size_t blockDelta = blockPos - it->first;
        return SeqPos(last.getSeqIndex(), last.getSeqPos() + blockDelta);
}

size_t SeqBlock::getRemainingSeqLen(size_t blockPos) const
{
        // check the validity of the input parameter
        assert (blockPos < block.size());

        // find the last marker less than or equal to blockPos
        auto it = block2seq.upper_bound(blockPos);
        assert(it != block2seq.begin());

        size_t nxtBlockPos = (it == block2seq.end()) ? block.size() : it->first;
        return nxtBlockPos - blockPos;
}

void SeqBlock::append(const std::string& data, const SeqPos& sp)
{
        // if necessary, insert a new marker
        if (!isContiguous(sp))
                block2seq[block.size()] = sp;

        block.append(data);
}

void SeqBlock::getSuffixBlock(SeqBlock& sb, size_t startPos)
{
        sb.clear();
        sb.block = block.substr(startPos);

        // find the last marker less than or equal to startPos
        auto it = block2seq.upper_bound(startPos);
        assert(it != block2seq.begin());
        --it;

        // copy the primary sequence marker
        size_t blockDelta = startPos - it->first;
        sb.block2seq[0] = SeqPos(it->second.getSeqIndex(),
                                 it->second.getSeqPos() + blockDelta);

        // copy the remainder of sequence markers (if any)
        for (it++; it != block2seq.end(); it++)
                sb.block2seq[it->first - startPos] = it->second;
}

std::ostream& operator<< (std::ostream& os, const SeqBlock& sb)
{
        os << sb.block << "\n";
        for (auto el : sb.block2seq)
                os << el.first << ": (" << el.second.getSeqIndex()
                   << ", " << el.second.getSeqPos() << ")" << "\n";
        return os;
}

// ============================================================================
// FASTA BATCH HANDLER
// ============================================================================

bool FastaBatch::moveToNextFile()
{
        // close input stream if necessary
        if (_ifs.is_open())
                _ifs.close();

        // return false if no files are left
        if (_currFileIdx >= seqFiles.size())
                return false;

        // open new sequence file
        _ifs.open(seqFiles[_currFileIdx]);
        if (!_ifs)
                throw runtime_error("Could not open file: " + seqFiles[_currFileIdx]);

        _currFileIdx++;
        return true;
}

bool FastaBatch::getNextLine(std::string& line, SeqPos& seqPos)
{
        while (true) {
                // move to the next file if necessary
                if (!getline(_ifs, line))
                        if (!moveToNextFile())
                                return false;

                // skip empty lines
                if (line.empty())
                        continue;

                // skip fasta sequence descriptor lines
                if (line.front() == '>') {
                        _currSeqLen = 0;
                        istringstream iss(line.substr(1));
                        string seqName;
                        iss >> seqName;
                        _seqNames.push_back(seqName);
                        continue;
                }

                // at this point we have a non-empty sequence line
                break;
        }

        // check for fasta format
        if (_seqNames.empty())
                throw runtime_error("Input file does not appear to be in fasta format\n");

        // set the sequence position
        seqPos = SeqPos(_seqNames.size() - 1, _currSeqLen);

        _totSeqLen += line.size();
        _currSeqLen += line.size();

        return true;
}

void FastaBatch::filterLine(const string& line, const SeqPos& seqPos,
                            deque<std::string>& filtLine,
                            deque<SeqPos>& filtSeqPos) const
{
        SeqPos currSeqPos = seqPos;
        bool isOpen = false;
        size_t openIdx = 0;

        for (size_t i = 0; i < line.size(); i++) {
                bool isValid = ((line[i] == 'A') || (line[i] == 'a') ||
                                (line[i] == 'C') || (line[i] == 'c') ||
                                (line[i] == 'G') || (line[i] == 'g') ||
                                (line[i] == 'T') || (line[i] == 't'));

                // open a sequence stream
                if (isValid && !isOpen) {
                        filtSeqPos.push_back(currSeqPos);
                        isOpen = true;
                        openIdx = i;
                }

                // close a sequence stream
                if (!isValid && isOpen) {
                        filtLine.push_back(line.substr(openIdx, i - openIdx));
                        isOpen = false;
                }

                currSeqPos.incSeqPos();
        }

        // also store the trailing line
        if (isOpen)
                filtLine.push_back(line.substr(openIdx, line.size() - openIdx));
}

bool FastaBatch::appendNextBlock(SeqBlock& block, size_t maxSize)
{
        bool retVal = false;
        // reduce maxSize in order not to exceed to total maximum filtered content size
        maxSize = min<size_t>(maxSize, _maxFiltSeqLen - _totFiltSeqLen + block.size());

        while (block.size() < maxSize) {
                // get more sequence data if necessary
                while (_filtSeqv.empty()) {
                        string tmp; SeqPos tmpSeqPos;
                        if (!getNextLine(tmp, tmpSeqPos))
                                return retVal;
                        filterLine(tmp, tmpSeqPos, _filtSeqv, _filtSeqPosv);
                }

                // copy filtered sequence data onto seqPos
                size_t thisSize = min<size_t>(_filtSeqv.front().size(),
                                              maxSize - block.size());
                block.append(_filtSeqv.front().substr(0, thisSize),
                             _filtSeqPosv.front());
                _totFiltSeqLen += thisSize;
                retVal = true;

                // update or remove the _filtSeq data
                if (_filtSeqv.front().size() == thisSize) {
                        _filtSeqv.pop_front();
                        _filtSeqPosv.pop_front();
                } else {
                        _filtSeqv.front() = _filtSeqv.front().substr(thisSize);
                        _filtSeqPosv.front().incSeqPos(thisSize);
                }
        }

        return retVal;
}

bool FastaBatch::getNextFilteredLine(std::string& line, SeqPos& seqPos)
{
        // this function is thread-safe
        lock_guard<mutex> lock(m);

        // get more sequence data if necessary
        while (_filtSeqv.empty()) {
                string tmp; SeqPos tmpSeqPos;
                if (!getNextLine(tmp, tmpSeqPos))
                        return false;
                filterLine(tmp, tmpSeqPos, _filtSeqv, _filtSeqPosv);
        }

        line = _filtSeqv.front();
        seqPos = _filtSeqPosv.front();
        _totFiltSeqLen += line.size();
        _filtSeqv.pop_front();
        _filtSeqPosv.pop_front();

        return true;
}

bool FastaBatch::getNextOverlappingBlock(SeqBlock& block, size_t payload,
                                         size_t overlap)
{
        // this function is thread-safe
        lock_guard<mutex> lock(m);

        // first copy the overlap from the previous block
        block = _nextBlock;

        // reduce maxSize in order not to exceed the total maximum size
        appendNextBlock(block, payload + overlap);

        // copy the overlap to _nextBlock to include as payload in future calls
        if (block.size() > payload)
                block.getSuffixBlock(_nextBlock, payload);
        else
                _nextBlock.clear();

        return !block.empty();
}

// ============================================================================
// SEQUENCE MATRIX
// ============================================================================

bool SeqMatrix::getNextSeqMatrix(FastaBatch& bf)
{
        // get the next block
        if(!bf.getNextOverlappingBlock(block, K * numRow, overlap))
                return false;

        // fill the sequence matrix
        S.setZero();
        for (size_t i = 0; i < numRow; i++) {
                // we're in the payload section of S
                for (size_t j = 0; j < K; j++) {
                        if ((i*K+j) >= block.size())
                                return true;
                        if (block[i*K+j] == 'A')
                                S(i, 4*j+0) = 1;
                        if (block[i*K+j] == 'C')
                                S(i, 4*j+1) = 1;
                        if (block[i*K+j] == 'G')
                                S(i, 4*j+2) = 1;
                        if (block[i*K+j] == 'T')
                                S(i, 4*j+3) = 1;
                        numOccRow = i + 1;
                }

                // we're in the overlap section of S -- don't return (!)
                for (size_t j = K; j < (K + overlap); j++) {
                        if ((i*K+j) >= block.size())
                                break;
                        if (block[i*K+j] == 'A')
                                S(i, 4*j+0) = 1;
                        if (block[i*K+j] == 'C')
                                S(i, 4*j+1) = 1;
                        if (block[i*K+j] == 'G')
                                S(i, 4*j+2) = 1;
                        if (block[i*K+j] == 'T')
                                S(i, 4*j+3) = 1;
                        numOccRow = i + 1;
                }
        }

        return true;
}

std::ostream& operator<< (std::ostream& os, const SeqMatrix& sm)
{
        os << "Matrix dims: " << sm.S.nRows() << " x " << sm.S.nCols()
           << ", of which " << sm.numOccRow << " are occupied\n";
        os << sm.block << "\n";
        sm.S.printSequence(sm.overlap);
        return os;
}
