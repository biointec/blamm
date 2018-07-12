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

#ifndef ORTHO_H
#define ORTHO_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

// ============================================================================
// ORTHO COUNT CLASS
// ============================================================================

class OrthoCount
{
private:
        std::string motifName;
        std::vector<size_t> motifCounts;
        std::vector<std::vector<size_t> > randomCounts;

public:
        /**
         * Default constructor
         */
        OrthoCount(const std::string& motifName) : motifName(motifName) {}

        /**
         * Set the conserved counts for the true motif
         * @param counts Vector with conserved counts per BLS threshold
         */
        void setMotifCounts(const std::vector<size_t>& counts) {
                motifCounts = counts;
        }

        /**
         * Add the random counts for the true motif
         * @param counts Vector with conserved counts per BLS threshold
         */
        void addRandomCounts(const std::vector<size_t>& counts) {
                randomCounts.push_back(counts);
        }

        size_t size() const {
                return randomCounts.size();
        }

        std::string getName() const {
                return motifName;
        }

        std::vector<float> getRandomMedian() const {
                std::vector<float> median;
                for (size_t i = 0; i < motifCounts.size(); i++) {

                        std::vector<size_t> randoms;
                        for (size_t j = 0; j < randomCounts.size(); j++)
                                randoms.push_back(randomCounts[j][i]);

                        std::sort(randoms.begin(), randoms.end());

                        float randomCount;
                        if (randoms.size() == 0) {
                                randomCount = 0;
                        } else if (randoms.size() % 2 == 1) {
                                randomCount = randoms[randoms.size()/2];
                        } else {   // size is even
                                float el1 = randoms[randoms.size()/2 - 1];
                                float el2 = randoms[randoms.size()/2];
                                randomCount = 0.5 * (el1 + el2);
                        }

                        median.push_back(randomCount);
                }

                return median;
        }

        std::vector<float> computeCScores() {

                std::vector<float> CScores;
                std::vector<float> medianCounts = getRandomMedian();

                for (size_t i = 0; i < motifCounts.size(); i++) {
                        float C = 1.0 - medianCounts[i] / motifCounts[i];
                        if (motifCounts[i] == 0)
                                C = 0;
                        if (C < 0)
                                C = 0;

                        CScores.push_back(C);
                }

                return CScores;
        }

        void print();
};

// ============================================================================
// ORTHOGROUP
// ============================================================================

class OrthoGroup
{
private:
        std::set<std::string> sequences;
        std::set<std::string> species;

public:
        OrthoGroup() {}

        void insert(const std::string& speciesName, const std::string& geneName) {
                sequences.insert(geneName);
                species.insert(speciesName);
        }

        size_t size() const {
                return sequences.size();
        }

        size_t getNumSpecies() const {
                return species.size();
        }

        std::set<std::string>::const_iterator seqBegin() const {
                return sequences.begin();
        }

        std::set<std::string>::const_iterator seqEnd() const {
                return sequences.end();
        }

        const std::set<std::string>& getSpecies() const {
                return species;
        }

        /**
         * Print ortho group to the screen
         * @param os Output stream to add to
         * @param og Ortho group to print
         * @return Output stream with the matrix elements
         */
        friend std::ostream& operator<< (std::ostream& os, const OrthoGroup& og);
};

// ============================================================================
// ORTHOCONTAINER
// ============================================================================

class OrthoContainer
{
private:
        std::map<std::string, OrthoGroup> orthoGroups;
        std::multimap<std::string, std::string> seq2ortho;

        typedef std::multimap<std::string, std::string>::iterator OrthIt;

public:
        /**
         * Default Constructor
         */
        OrthoContainer() {}

        void makeHistogram(size_t numSpecies);

        /**
         * Load the ortho container from disk
         * @param filename filename of the input file
         */
        void load(const std::string& filename);

        /**
         * Return the size of the container
         * @return The size of the container
         */
        size_t size() const {
                return orthoGroups.size();
        }

        const OrthoGroup& getOrthoGroup(const std::string & orthoName) const {
                return orthoGroups.at(orthoName);
        }

        /**
         * Get a range of ortho groups that contain a specific sequence
         * @return Range of ortho groups
         */
        std::pair<OrthIt, OrthIt> equal_range(std::string& seqName) {
                return seq2ortho.equal_range(seqName);
        }
};

// ============================================================================
// ORTHO MODULE
// ============================================================================

class Ortho
{
private:
        /**
         * Print module instructions
         */
        void printUsage() const;

public:
        /**
         * Constructor (run Ortho module)
         * @param argc Command line argument count
         * @param argv Command line argument values
         */
        Ortho(int argc, char **argv);
};

#endif
