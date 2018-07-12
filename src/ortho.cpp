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

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "ortho.h"
#include "phylotree.h"
#include "species.h"
#include "motif.h"

using namespace std;

void OrthoCount::print()
{
        cout << motifName << endl;
        for (auto it : motifCounts)
                cout << it << " ";
        cout << "\n";
        vector<float> median = getRandomMedian();
        for (auto it : median)
                cout << it << " ";
        cout << "\n";
        vector<float> C = computeCScores();
        for (auto it : C)
                cout << it << " ";
        cout << "\n";
}

// ============================================================================
// ORTHOCONTAINER
// ============================================================================

void OrthoContainer::load(const string& filename)
{
        // read the input file and create the orthoGroups
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        while (true) {
                string refSpecies, refSeqName, orthoSpecies, orthoSeqName;
                ifs >> refSpecies >> refSeqName >> orthoSpecies >> orthoSeqName;
                if (!ifs)
                        break;

                orthoGroups[refSeqName].insert(refSpecies, refSeqName);
                orthoGroups[refSeqName].insert(orthoSpecies, orthoSeqName);
        }

        // create a sequence 2 ortho index
        for (auto it = orthoGroups.begin(); it != orthoGroups.end(); it++) {
                const string& orthoName = it->first;
                const OrthoGroup& og = it->second;
                for (auto it2 = og.seqBegin(); it2 != og.seqEnd(); it2++)
                        seq2ortho.insert(make_pair(*it2, orthoName));
        }
}

void OrthoContainer::makeHistogram(size_t numSpecies)
{
        vector<size_t> hist(numSpecies+1, 0);

        for (const auto& it : orthoGroups) {
                const OrthoGroup& og = it.second;
                hist[og.getNumSpecies()]++;
        }

        for (size_t i = 1; i <= numSpecies; i++)
                cout << "#orthogroups with " << i << " species: " << hist[i] << endl;

        /*vector<size_t> hist2(numSpecies+10, 0);

        for (const auto& it : orthoGroups) {
                const OrthoGroup& og = it.second;
                hist2[og.size()]++;
        }

        for (size_t i = 1; i <= numSpecies+10; i++)
                cout << "#orthogroups with " << i << " sequences: " << hist2[i] << endl;*/
}

// ============================================================================
// ORTHO MODULE
// ============================================================================

std::ostream& operator<< (std::ostream& os, const OrthoGroup& og)
{
        for (auto it : og.sequences)
                os << it << "\n";
        return os;
}

void Ortho::printUsage() const
{
        cout << "Usage: blstools ortho [options] motifs.input sequences.input orthogroups.input phylotree.input occurrences.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tfilename for output BLS values [default = stdout]\n\n";

        cout << " File \"phylotree.input\" should contain the phylogenetic tree in Newick format, e.g.:\n";
        cout << "  \"(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);\"\n\n";

        cout << " File \"leafs.input\" should contain a list of space or tab separated\n";
        cout << "  leaf names for which the BLS should be computed (one line per entry)\n\n";

        cout << " Example:\n";
        cout << "  blstools bls -o output.bls tree.newick leafs.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

Ortho::Ortho(int argc, char ** argv)
{
        // check for sufficient arguments
        if (argc < 7) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
        for (int i = 2; i < argc-5; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        string motifFilename = argv[argc-5];
        string manifestFilename = argv[argc-4];
        string orthoFilename = argv[argc-3];
        string phyloFilename = argv[argc-2];
        string occFilename = argv[argc-1];

        // A) Load the motifs
        MotifContainer motifContainer;
        motifContainer.load(motifFilename, true);
        cout << "Loaded " << motifContainer.size() << " motifs from disk\n";

        // B) Load the manifest file
        string dictFilename = manifestFilename + ".dict";

        SpeciesContainer speciesContainer;
        speciesContainer.load(dictFilename);
        cout << "Loaded dictionaire with " << speciesContainer.size() << " species\n";

        // C) Read the ortho groups
        OrthoContainer orthoContainer;
        orthoContainer.load(orthoFilename);
        cout << "Loaded " << orthoContainer.size() << " orthology groups\n";

        orthoContainer.makeHistogram(speciesContainer.size());

         // D) process the phylogenetic tree
        ifstream ifs(phyloFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + phyloFilename);
        string newickStr;
        getline(ifs, newickStr);
        ifs.close();

        PhylogeneticTree pt(newickStr);
        pt.normalizeBranchLength();

        // check whether the species in the phylotree match with the ones in the manifest file
        set<string> speciesNames = pt.getAllNames();
        for (auto it : speciesContainer)
                if (speciesNames.find(it.getName()) == speciesNames.end())
                        throw runtime_error("ERROR: Species name \"" + it.getName() + "\" does not occur in " + phyloFilename);

        cout << "Loaded phylogenetic tree" << endl;

        // E) Read the occurrence file
        ofstream ofsBLS("motifBLS.txt");
        ofstream ofsCounts("motifCounts.txt");

        ifs.open(occFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + occFilename);

        size_t numBLSIntv = 10;

        for (size_t i = 0; i < motifContainer.size(); i++) {
                map<string, set<string> > orthoSpecComb;
                map<string, set<string> > orthoGeneComb;

                while (ifs) {
                        // read a single line of the occurrence file
                        string temp;
                        size_t motifID, speciesID, seqID;
                        int currFP = ifs.tellg();
                        ifs >> motifID >> speciesID >> seqID >> temp >> temp >> temp;

                        if (!ifs)
                                break;          // eof reached

                        if (speciesID > speciesContainer.size())
                                throw runtime_error("ERROR: File " + occFilename + " contains a speciesID not present in file " + manifestFilename);

                        if (motifID > motifContainer.size())
                                throw runtime_error("ERROR: File " + occFilename + " contains a motifID not present in file " + motifFilename);

                        if (motifID < i)
                                throw runtime_error("ERROR: File " + occFilename + " is not sorted. Sort this file prior to running the ortho module");

                        // if the line deals with another motif, rewind
                        if (motifID != i) {
                                ifs.seekg(currFP, ios_base::beg);
                                break;
                        }

                        // it deals with the current motif, so score it
                        string species = speciesContainer.getSpecies(speciesID).getName();
                        string seq = speciesContainer.getSpecies(speciesID).getSeqName(seqID);

                        // for all ortho groups that contain sequence "seq"
                        auto range = orthoContainer.equal_range(seq);
                        for (auto it = range.first; it != range.second; it++) {
                                orthoSpecComb[it->second].insert(species);
                                orthoGeneComb[it->second].insert(seq);
                        }
                }

                // count the number of conserved instances for each BLS threshold
                vector<size_t> counts(numBLSIntv, 0);
                for (auto it : orthoSpecComb) {
                        float BLS = pt.getBLS(it.second);

                        // const OrthoGroup& og = orthoContainer.getOrthoGroup(it.first);
                        // float maxBLS = pt.getBLS(og.getSpecies());
                        // if (maxBLS > 0)
                        //       BLS /= maxBLS;

                        ofsBLS << i << "\t" << motifContainer[i].getName() << "\t" << it.first << "\t" << BLS;

                        for (auto gene : orthoGeneComb[it.first])
                                ofsBLS << "\t" << gene;
                        ofsBLS << endl;

                        float intWidth = 1.0 / numBLSIntv;
                        size_t end = floor(BLS / intWidth) + 1;
                        if (end > numBLSIntv)
                                end = numBLSIntv;

                        for (size_t i = 0; i < end; i++)
                                counts[i]++;
                }

                cout << i << "\t" << motifContainer[i].getName();
                ofsCounts << i << "\t" << motifContainer[i].getName();
                for (size_t j = 0; j < numBLSIntv; j++) {
                        cout << "\t" << counts[j];
                        ofsCounts << "\t" << counts[j];
                }
                cout << "\n";
                ofsCounts << "\n";
        }

        ifs.close();
        ofsBLS.close();
        ofsCounts.close();

        ofstream ofsCScores("motifCScores.txt");
        ifs.open("motifCounts.txt");
        if (!ifs)
                throw runtime_error("Could not open file: counts.txt");

        OrthoCount orthoCount("");
        while (ifs) {
                size_t motifID;
                vector<size_t> counts(numBLSIntv, 0);
                string motifName;

                ifs >> motifID >> motifName;
                for (size_t i = 0; i < numBLSIntv; i++)
                        ifs >> counts[i];

                // if we encounter a new motif
                if ((!ifs) || (!motifContainer[motifID].isPermutation())) {
                        if (orthoCount.size() > 0) {
                                ofsCScores << orthoCount.getName();
                                vector<float> C = orthoCount.computeCScores();
                                for (auto it : C)
                                        ofsCScores << "\t" << it;
                                ofsCScores << "\n";
                        }

                        orthoCount = OrthoCount(motifContainer[motifID].getName());
                        orthoCount.setMotifCounts(counts);
                } else {        // it's a random motif
                        orthoCount.addRandomCounts(counts);
                }
        }

        ofsCScores.close();
}
