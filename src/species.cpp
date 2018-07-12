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

#include <numeric>

#include "species.h"
#include "sequence.h"

using namespace std;

// ============================================================================
// SPECIES
// ============================================================================

void Species::computeNucleotideComposition()
{
        const size_t blockSize = 65536;

        int char2Idx[256];
        for (size_t i = 0; i < 256; i++)
                char2Idx[i] = 4;
        char2Idx[(unsigned short)'A'] = 0;
        char2Idx[(unsigned short)'a'] = 0;
        char2Idx[(unsigned short)'C'] = 1;
        char2Idx[(unsigned short)'c'] = 1;
        char2Idx[(unsigned short)'G'] = 2;
        char2Idx[(unsigned short)'g'] = 2;
        char2Idx[(unsigned short)'T'] = 3;
        char2Idx[(unsigned short)'t'] = 3;

        nuclCounts.fill(0);

        vector<string> sequences(fastaFiles.begin(), fastaFiles.end());
        FastaBatch batch(sequences);
        SeqBlock block;

        while (batch.getNextOverlappingBlock(block, blockSize, 0))
                for (size_t i = 0; i < block.size(); i++)
                        nuclCounts[char2Idx[(unsigned short)block[i]]]++;

        totSeqLen = accumulate(nuclCounts.begin(), nuclCounts.end(), 0ull);

        seqNames = batch.getSeqNames();
}

std::array<float, 4> Species::getNuclProbabilities(float pseudoCount) const
{
        float bgTotCounts = (float)accumulate(nuclCounts.begin(), nuclCounts.end(), 0ull);
        bgTotCounts += 4.0f * pseudoCount;
        std::array<float, 4> bgProb;
        for (size_t i = 0; i < 4; i++)
                bgProb[i] = (float(nuclCounts[i]) + pseudoCount) / bgTotCounts;
        return bgProb;
}

void Species::writeNuclProbabilities(float pseudoCount) const
{
        array<float, 4> bgProb = getNuclProbabilities(pseudoCount);
        auto oldPrec = cout.precision();
        cout.precision(3);
        cout << " [A: " << 100.0f*bgProb[0] << "%, C: " << 100.0f*bgProb[1]
             << "%, G: " << 100.0f*bgProb[2] << "%, T: " << 100.0f*bgProb[3]
             << "%]\n";
        cout.precision(oldPrec);
}

void Species::write(std::ofstream& ofs) const
{
        ofs << "SPECIES" << "\t" << name << "\n";
        ofs << "NUM_FASTA_FILES" << "\t" << fastaFiles.size() << "\n";

        for (const auto& file : fastaFiles)
                ofs << file << "\n";

        ofs << "TOT_SEQ_LENGTH" << "\t" << getTotalSeqLength() << "\n";

        ofs << "NUCL_COUNT_ACGT" << "\t"
            << nuclCounts[0] << "\t"
            << nuclCounts[1] << "\t"
            << nuclCounts[2] << "\t"
            << nuclCounts[3] << "\n";

        ofs << "NUM_SEQUENCES" << "\t" << seqNames.size() << "\n";
        for (auto& seqName : seqNames)
                ofs << seqName << "\n";
}

void Species::load(std::ifstream& ifs)
{
        string temp;

        // species name
        ifs >> temp >> name;
        assert(temp == "SPECIES");

        // number of fasta files
        size_t numFastaFiles;
        ifs >> temp >> numFastaFiles;
        assert(temp == "NUM_FASTA_FILES");

        // fasta files
        for (size_t j = 0; j < numFastaFiles; j++) {
                ifs >> temp;
                fastaFiles.insert(temp);
        }

        // total sequence length
        ifs >> temp >> totSeqLen;
        assert(temp == "TOT_SEQ_LENGTH");

        // nucleotide composition
        ifs >> temp >> nuclCounts[0] >> nuclCounts[1] >> nuclCounts[2]
            >> nuclCounts[3];
        assert(temp == "NUCL_COUNT_ACGT");

        // sequence names
        size_t numSequences;
        ifs >> temp >> numSequences;
        assert(temp == "NUM_SEQUENCES");
        seqNames.clear();
        seqNames.reserve(numSequences);
        for (size_t i = 0; i < numSequences; i++) {
                ifs >> temp;
                seqNames.push_back(temp);
        }
}

// ============================================================================
// SPECIES CONTAINER
// ============================================================================

void SpeciesContainer::addSeqFile(const string& speciesName,
                                  const string& fastaFilename)
{
        // create a new species if necessary
        size_t idx;
        auto it = name2Idx.find(speciesName);
        if (it == name2Idx.end()) {
                idx = speciesContainer.size();
                name2Idx[speciesName] = idx;
                speciesContainer.push_back(Species(speciesName));
        } else {
                idx = it->second;
        }

        speciesContainer[idx].addFile(fastaFilename);
}

void SpeciesContainer::computeNucleotideComposition()
{
        for (auto& species : speciesContainer) {
                cout << "Computing nucleotide composition for " << species.getName() << " ...";
                cout.flush();
                species.computeNucleotideComposition();
                cout << "\n";
        }
}

void SpeciesContainer::write(const std::string& filename)
{
        ofstream ofs(filename.c_str());
        if (!ofs)
                throw runtime_error("Cannot write to file: " + filename);

        ofs << "NUM_SPECIES" << "\t" << speciesContainer.size() << "\n";

        for (auto it : speciesContainer)
                it.write(ofs);

        ofs.close();
}

void SpeciesContainer::load(const std::string& filename)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Cannot open file: " + filename);

        // number of species
        string temp; size_t numSpecies;
        ifs >> temp >> numSpecies;
        assert(temp == "NUM_SPECIES");

        for (size_t i = 0; i < numSpecies; i++) {
                speciesContainer.push_back(Species("dummy"));
                speciesContainer.back().load(ifs);
        }

        ifs.close();
}
