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

#ifndef PWM_H
#define PWM_H

#include <vector>
#include "settings.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Species;
class MotifContainer;
class FastaBatch;
class ScoreHistogram;
class SeqMatrix;
class Matrix;

// ============================================================================
// HISTOGRAM MODULE
// ============================================================================

class Histogram
{
private:
        /**
         * Print module instructions
         */
        void printUsage() const;

        void extractObsScore(const Matrix& R, size_t offset,
                             const SeqMatrix& sm,
                             const MotifContainer& motifContainer,
                             std::vector<ScoreHistogram>& histContainer);

        void histThread(const MotifContainer& motifs,
                        FastaBatch& fb, std::vector<ScoreHistogram>& histContainer);

        void generateEmpiricalHist(const Species& species,
                                   const MotifContainer& MotifContainer,
                                   std::vector<ScoreHistogram>& histContainer);

        void generateTheoreticalHist(const Species& species,
                                     const MotifContainer& MotifContainer,
                                     std::vector<ScoreHistogram>& histContainer);

        Settings settings;      // settings object
        size_t maxLength;       // maximum length of sequence data to analyze
        size_t numBins;         // number of bins per histogram
        std::string histdir;    // histogram directory
        size_t numThreads;      // number of threads
        bool empirical;         // compute empirical spectra

public:
        /**
         * Constructor (run Histogram module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        Histogram(int argc, char **argv);
};

#endif
