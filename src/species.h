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

#ifndef SPECIES_H
#define SPECIES_H

#include <string>
#include <set>
#include <vector>
#include <map>

// ============================================================================
// SPECIES
// ============================================================================

class Species {
private:
        std::string name;                       // name of the species
        std::set<std::string> fastaFiles;       // fastaFiles

        std::array<size_t, 4> nuclCounts;       // nucleotide counts ACGT
        size_t totSeqLen;                       // total sequence length

        std::vector<std::string> seqNames;      // sequence names

public:
        /**
         * Default constructor
         */
        Species(std::string name) : name(name) {};

        /**
         * Get the name of the species
         * @return The name of the species
         */
        std::string getName() const {
                return name;
        }

        /**
         * Get the name of a sequence
         * @param seqID Sequence identifier
         * @return The name of the corresponding sequence
         */
        std::string getSeqName(size_t seqID) const {
                return seqNames[seqID];
        }

        /**
         * Compute the background frequencies for this species
         */
        void computeNucleotideComposition();

        /**
         * Add a sequence fasta file to this species
         * @param fastaFilename Fasta filename
         */
        void addFile(const std::string& fastaFilename) {
                fastaFiles.insert(fastaFilename);
        }

        /**
         * Get the sequence filename
         * @return A vector containing the sequence filenames
         */
        std::vector<std::string> getSequenceFilenames() const {
                return std::vector<std::string>(fastaFiles.begin(),
                                                fastaFiles.end());
        }

        /**
         * Get the nucleotide counts
         * @return the nucleotide counts
         */
        std::array<size_t, 4> getNuclCounts() const {
                return nuclCounts;
        }

        /**
         * Get the nucleotide probabilities
         * @param pseudoCount pseudo count to include when computing probabilities
         * @return the nucleotide probabilities
         */
        std::array<float, 4> getNuclProbabilities(float pseudoCount) const;

        /**
         * Write the nucleotide probabilities to stdout
         * @param pseudoCount pseudo count to include when computing probabilities
         */
        void writeNuclProbabilities(float pseudoCount) const;

        /**
         * Get the total sequence length
         * @return The total sequence length
         */
        size_t getTotalSeqLength() const {
                return totSeqLen;
        }

        /**
         * Write object to disk
         * @param ofs Opened and valid output filestream
         */
        void write(std::ofstream& ofs) const;

        /**
         * Load object from disk
         * @param ifs Opened and valid input filestream
         */
        void load(std::ifstream& ifs);
};

// ============================================================================
// SPECIES CONTAINER
// ============================================================================

class SpeciesContainer {
private:
        std::vector<Species> speciesContainer;
        std::map<std::string, size_t> name2Idx;

public:
        /**
         * Default constructor
         */
        SpeciesContainer() {};

        /**
         * Add a < species name, sequence file > to the container
         * @param speciesName Species name
         * @param fastaFilename Fasta filename
         */
        void addSeqFile(const std::string& speciesName,
                        const std::string& fastaFilename);

        /**
         * Compute the background frequencies for all species
         */
        void computeNucleotideComposition();

        /**
         * Write object to disk
         * @param filename Dictionary output filename
         */
        void write(const std::string& filename);

        /**
         * Load object from disk
         * @param filename Dictionary input filename
         */
        void load(const std::string& filename);

        /**
         * Get a species by identifier
         * @return A const-ref to the species
         */
        const Species& getSpecies(size_t speciesID) const {
                return speciesContainer[speciesID];
        }

        /**
         * Operator [] overloading
         * @return A const-ref to the species
         */
        const Species& operator[](size_t speciesID) const {
                return speciesContainer[speciesID];
        }

        /**
         * Get a species by name
         * @return A const-ref to the species
         */
        const Species& getSpecies(const std::string& speciesName) const {
                return speciesContainer[name2Idx.at(speciesName)];
        }

        /**
         * Return an iterator pointing to the first species in the container
         * @return An iterator to the beginning of the sequence container
         */
        std::vector<Species>::const_iterator begin() const {
                return speciesContainer.begin();
        }

        /**
         * Return an iterator pointing past the final species in the container
         * @return An iterator to the end of the sequence container
         */
        std::vector<Species>::const_iterator end() const {
                return speciesContainer.end();
        }

        /**
         * Get the number of species
         * @return The number of species
         */
        size_t size() const {
                return speciesContainer.size();
        }
};

#endif
