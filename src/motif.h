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

#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <array>
#include <atomic>
#include <cassert>
#include <map>

#include "matrix.h"

// ============================================================================
// SCORE HISTOGRAM
// ============================================================================

class ScoreHistogram
{
private:
        std::atomic<size_t>* counts;          // score counts
        float minScore;         // minimum score of the histogram
        float maxScore;         // maximum score of the histogram
        size_t numBins;         // number of bins
        float width;            // width of a bin

public:
        /**
         * Default constructor
         */
        ScoreHistogram() : counts(NULL), minScore(0.0f), maxScore(0.0f),
                numBins(0), width(0.0f) {}

        /**
         * Default constructor
         * @param minScore Minimum score for the histogram
         * @param maxScore Maximum score for the histogram
         * @param numBins Number of bins in the histogram
         */
        ScoreHistogram(float minScore, float maxScore, size_t numBins) :
                minScore(minScore), maxScore(maxScore), numBins(numBins)
        {
                assert(numBins > 0);

                width = (maxScore - minScore) / (float)numBins;
                counts = new std::atomic<size_t>[numBins];
                for (size_t i = 0; i < numBins; i++)
                        counts[i] = 0;
        }

        /**
         * Copy constructor
         * @param S Histogram to copy
         */
        ScoreHistogram(const ScoreHistogram& S) {
                minScore = S.minScore;
                maxScore = S.maxScore;
                numBins = S.numBins;
                width = S.width;

                counts = new std::atomic<size_t>[numBins];
                for (size_t i = 0; i < numBins; i++)
                        counts[i].store(S.counts[i].load());
        }

        /**
         * Destructor
         */
        ~ScoreHistogram() {
                if (counts != NULL)
                        delete [] counts;
        }

        /**
         * Add a single observation to the histogram
         * @param score Score of the observation
         */
        void addObservation(float score) {
                int binIdx = int((score - minScore) / width);
                binIdx = std::max<int>(0, binIdx);
                binIdx = std::min<int>(numBins-1, binIdx);

                counts[binIdx]++;
        }

        /**
         * Set the number of observations in the histogram
         * @param score Score of the observation
         * @param count Number of observations
         */
        void setNumObservations(float score, size_t count) {
                int binIdx = int((score - minScore) / width);
                binIdx = std::max<int>(0, binIdx);
                binIdx = std::min<int>(numBins-1, binIdx);

                counts[binIdx].store(count);
        }

        /**
         * Compute the weighted average of the observations
         */
        float getAverage() const;

        /**
         * Compute the score cutoff corresponding to a certain pvalue
         * @param pvalue p-value
         * @return The score cutoff
         */
        float getScoreCutoff(float pvalue) const;

        /**
         * Write a GNUplot file containing the histogram
         * @param dir Output directory (must end by '/')
         * @param baseFilename Base filename (.dat and .gnu will be added)
         * @param label Histogram label string
         */
        void writeGNUPlotFile(const std::string& dir,
                              const std::string& baseFilename,
                              const std::string& label) const;

        /**
         * Load a histogram from disk
         * @param dir Output directory (must end by '/')
         * @param baseFilename Base filename (.dat will be added)
         */
        void loadHistogram(const std::string& dir,
                           const std::string& baseFilename);
};

// ============================================================================
// MOTIF
// ============================================================================

class Motif {
private:
        std::string name;                               // name of the motif
        std::vector<std::array<size_t, 4> > PFM;        // position frequency matrix
        std::vector<std::array<float, 4> > PPM;         // position probability matrix
        std::vector<std::array<float, 4> > PWM;         // position weight matrix
        float threshold;                                // occurrence threshold cutoff
        bool revComp;                                   // is this a rev-compl motif

public:
        void computeTheoreticalSpectrum(size_t numBins,
                                        const std::array<float, 4>& background,
                                        std::map<float, float>& spectrum) const;

        /**
         * Default constructor
         */
        Motif() : threshold(0.0f), revComp(false) {}

        /**
         * Default constructor
         * @param name Name of the motif
         */
        Motif(const std::string& name) : name(name),
                threshold(0.0f), revComp(false) {}

        /**
         * Get the name of the motif
         * @return The name of the motif
         */
        const std::string& getName() const {
                return name;
        }

        /**
         * Is the motif a reverse complement?
         * @return True or false
         */
        bool isRevCompl() const {
                return revComp;
        }

        /**
         * Get the motif score threshold
         * @return Motif score threshold
         */
        float getThreshold() const {
                return threshold;
        }

        /**
         * Add a character to the positoin frequency matrix (PFM)
         * @param c counts of the character to add
         */
        void addCharacter(const std::array<size_t, 4>& c) {
                PFM.push_back(c);
        }

        /**
         * Convert a PFM into a PWM given background nucleotide counts
         * @param bgCounts Background nucleotide counts (ACTG)
         * @param pseudoCounts Pseudo counts for both the motif PFM and bgCounts
         */
        void PFM2PWM(const std::array<size_t, 4>& bgCounts,
                     float pseudoCounts);

        /**
         * Get the size (= length) of a motif
         * @return The size of a motif
         */
        size_t size() const {
                return PFM.size();
        }

        /**
         * Get the observation count in the PFM
         * @return The observation count
         */
        size_t getCount() const {
                if (PFM.empty())
                        return 0;
                return PFM[0][0] + PFM[0][1] + PFM[0][2] + PFM[0][3];
        }

        /**
         * Operator[] overloading to access the PWM
         * @param i position index
         */
        const std::array<float,4>& operator[](size_t i) const {
                return PWM[i];
        }

        /**
         * Get the position frequency matrix
         * @return the position frequency matrix
         */
        std::vector<std::array<size_t, 4> > getPFM() const {
                return PFM;
        }

        /**
         * Get the base name of the motif (with __XX for random permutations)
         * @return The base name of the motif
         */
        std::string getBaseName() const {
                std::string retval = name;
                ssize_t i = name.size()-1;
                for ( ; i >= 0; i--)
                        if (!isdigit(name[i]))
                                break;
                if (i < 1)
                        return retval;
                if ((name[i-1] != '_') || (name[i] != '_'))
                        return retval;
                return retval.substr(0, i-1);
        }

        /**
         * Check whether a motif is a permutation (__XX at the end)
         * @return true or false
         */
        bool isPermutation() const {
                ssize_t i = name.size()-1;
                for ( ; i >= 0; i--)
                        if (!isdigit(name[i]))
                                break;
                if (i < 1)
                        return false;
                return ((name[i-1] == '_') && (name[i] == '_'));
        }

        /**
         * Reverse complement the motif (both PFM and PWM)
         */
        void revCompl();

        /**
         * Write the MOODS file
         * @param filename Filename
         */
        void writeMOODSFile(const std::string& filename) const;

        /**
         * Get the maximum PWM score for this motif
         * @return The maximum PWM score for this motif
         */
        float getMaxScore() const;

        /**
         * Get the minimum PWM score for this motif
         * @return The minimum PWM score for this motif
         */
        float getMinScore() const;

        /**
         * Compute the PWM score given a pattern
         * @param pattern Pattern to compare against
         */
        float getScore(const std::string& pattern) const;

        /**
         * Set the motif score threshold
         * @param threshold_ Motif score threshold
         */
        void setThreshold(float threshold_) {
                threshold = threshold_;
        }

        /**
         * Operator< overloading based on size
         * @return true of false
         */
        bool operator<(const Motif& rhs) const {
                return size() < rhs.size();
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param m Motif
         */
        friend std::ostream& operator<< (std::ostream& os, const Motif& m);
};

// ============================================================================
// MATRIX CHUNK
// ============================================================================

class MatrixTile {
public:
        size_t rowStart;
        size_t colStart;
        size_t rowEnd;
        size_t colEnd;

        /**
         * Default constructor
         */
        MatrixTile() : rowStart(0), colStart(0), rowEnd(0), colEnd(0) {}

        /**
         * Constructor with element initialization
         * @param rowStart start row
         * @param colStart start column
         * @param rowEnd end row
         * @param colEnd end column
         */
        MatrixTile(size_t rowStart, size_t colStart, size_t rowEnd,
                   size_t colEnd) : rowStart(rowStart), colStart(colStart),
                   rowEnd(rowEnd), colEnd(colEnd) {}

        /**
         * Operator < overloading
         * @param rhs Right hand side object
         * @return True of false
         */
        bool operator< (const MatrixTile& rhs) {
                return rowStart < rhs.rowStart;
        }

        size_t getArea() const {
                return (rowEnd-rowStart) * (colEnd-colStart);
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param m Matrix tile
         */
        friend std::ostream& operator<< (std::ostream& os, const MatrixTile& m) {
                os << "motif [" << m.colStart << " - " << m.colEnd << "[ -- max. length = " << m.rowEnd;

                /*"(" << m.rowStart << "-" << m.colStart << ") -- "
                   << "(" << m.rowEnd << "-" << m.colEnd << ")";*/
                return os;
        }
};

// ============================================================================
// MOTIF CONTAINER
// ============================================================================

typedef std::pair<MatrixTile, MatrixTile> TilePair;

class MotifContainer {
private:
        /**
         * Given an input matrix tile, compute the optimal split
         * @param input Input matrix tile
         * @return Two output tiles
         */
        TilePair findBestSplit(const MatrixTile& input);

        /**
         * Keep/reject a split
         * @param tilePair Two input tiles
         * @param tileZeroMinArea Minimum zero area gained by splitting
         * @return True or false
         */
        bool keepSplit(const TilePair& tilePair, size_t tileMinZeroArea);

        std::vector<Motif> motifs;              // actual motifs
        Matrix P;                               // pattern matrix
        std::vector<MatrixTile> matrixTiles;    // block structor of P
        std::vector<size_t> col2MotifID;        // motif ID at column j in P

public:
        /**
         * Default constructor
         */
        MotifContainer() {}

        /**
         * Load a motif file from disk
         * @param filename Filename of the motif input file
         * @param loadPermutations Also load the random permutations
         */
        void load(const std::string& filename, bool loadPermutations);

        /**
         * Load cluster-buster motifs from disk
         * @param filename Filename of the motif input file
         * @param motifs Vector to store motifs (output)
         */
        void loadCBMotifs(const std::string& filename,
                          std::vector<Motif>& motifs);

        /**
         * Load Jaspar motifs from disk
         * @param filename Filename of the motif input file
         * @param motifs Vector to store motifs (output)
         */
        void loadJasparMotifs(const std::string& filename,
                              std::vector<Motif>& motifs);

        /**
         * Add the reverse complement motifs to the container
         */
        void addReverseComplements();

        /**
         * Fill in the elements of matrix P
         * @param bgCounts Background nucleotide counts (ACTG)
         * @param pseudoCount Pseudo counts for both the motif PFM and bgCounts
         */
        void generateMatrix(const std::array<size_t, 4>& bgCounts,
                            float pseudoCount);

        /**
         * Partition matrix P into tiles if needed
         * @param tileMinZeroArea Minimum area for a tile
         */
        void generateMatrixTiles(size_t tileMinZeroArea);

        /**
         * Get a const-reference to matrix P
         * @return A const-reference to matrix P
         */
        const Matrix& getMatrix() const {
                return P;
        }

        /**
         * Get a const-reference to the matrix tiles
         * @return A const-reference to the matrix tiles
         */
        const std::vector<MatrixTile>& getMatrixTiles() const {
                return matrixTiles;
        }

        /**
         * Get the motif ID contained in column i of the matrix
         * @return The motif ID at column i of the matrix
         */
        size_t getMotifIDAtCol(size_t i) const {
                return col2MotifID[i];
        }

        /**
         * Write the motif names
         * @param filename File name of the motif file
         */
        void writeMotifNames(const std::string& filename);

        /**
         * Write the possum file
         * @param filename File name of the possum file
         */
        void writePossumFile(const std::string& filename);

        /**
         * Write the MOODS file
         */
        void writeMOODSFiles();

        /**
         * Get the number of motifs
         */
        size_t size() const {
                return motifs.size();
        }

        const Motif& operator[](size_t index) const {
                return motifs[index];
        }

        size_t getMaxMotifLen() const;

        /**
         * Return an iterator pointing to the first motif in the container
         * @return An iterator to the beginning of the motif container
         */
        std::vector<Motif>::const_iterator begin() const {
                return motifs.begin();
        }

        /**
         * Return an iterator pointing past the final motif in the container
         * @return An iterator to the end of the motif container
         */
        std::vector<Motif>::const_iterator end() const {
                return motifs.end();
        }

        /**
         * Return an iterator pointing to the first motif in the container
         * @return An iterator to the beginning of the motif container
         */
        std::vector<Motif>::iterator begin() {
                return motifs.begin();
        }

        /**
         * Return an iterator pointing past the final motif in the container
         * @return An iterator to the end of the motif container
         */
        std::vector<Motif>::iterator end() {
                return motifs.end();
        }
};

// ============================================================================
// MOTIF OCCURRENCES
// ============================================================================

class MotifOccurrence {
private:
        unsigned int motifID;
        unsigned int speciesID;
        unsigned int sequenceID;
        size_t sequencePos;
        char strand;
        float score;

public:
        /**
         * Default constructor
         */
        MotifOccurrence() {};

        /**
         * Constructor with arguments
         * @param motifID The motif identifier
         * @param speciesID The species identifier
         * @param sequenceID The sequence identifier
         * @param sequencePos The sequence position identifier
         * @param strand + or - strand
         * @param score The motif score
         */
        MotifOccurrence(unsigned int motifID, unsigned int speciesID,
                        unsigned int sequenceID, size_t sequencePos,
                        char strand, float score) :
                motifID(motifID), speciesID(speciesID), sequenceID(sequenceID),
                sequencePos(sequencePos), strand(strand), score(score) {}

        /**
         * Get the motif identifier
         * @return The motif identifier
         */
        unsigned int getMotifID() const {
                return motifID;
        }

        /**
         * Get the species identifier
         * @return The species identifier
         */
        unsigned int getSpeciesID() const {
                return speciesID;
        }

        /**
         * Get the sequence identifier
         * @return The sequence identifier
         */
        unsigned int getSequenceID() const {
                return sequenceID;
        }

        /**
         * Get the sequence position
         * @return The sequence position
         */
        size_t getSequencePos() const {
                return sequencePos;
        }

        /**
         * Get the sequence strand
         * @return The sequence strand
         */
        char getStrand() const {
                return strand;
        }

        /**
         * Get the motif score
         * @return The motif score
         */
        float getScore() const {
                return score;
        }

        /**
         * operator<< overloading
         * @param os Output stream
         * @param m Motif occurrence
         */
        friend std::ostream& operator<< (std::ostream& os, const MotifOccurrence& m);
};

#endif
