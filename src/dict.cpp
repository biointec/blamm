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
#include <sstream>

#include "dict.h"
#include "sequence.h"
#include "species.h"

using namespace std;

void Dictionary::printUsage() const
{
        cout << "Usage: blamm dict [options] sequences.input\n";
        cout << "Goal: compute nucleotide frequencies for the input sequences "
             << "and generate a sequence dictionary file\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n\n";

        cout << " File \"sequences.input\" contains a list of input"
             << " fasta files in the following format:\n";
        cout << "   speciesID_1\tspecies1_sequences.fasta\n";
        cout << "   speciesID_2\tspecies2_sequences.fasta\n";
        cout << "   speciesID_3\tspecies2_sequences.fasta\n";
        cout << "   ...\n";
        cout << " where speciesID_x is a user-defined identifier per species\n\n";

        cout << " Example:\n";
        cout << "  blamm dict sequences.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

Dictionary::Dictionary(int argc, char ** argv)
{
        // check for sufficient arguments
        if (argc < 3) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
        for (int i = 2; i < argc-1; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        cout << "Welcome to blamm -- dictionary model" << endl;

        // read the manifest file
        SpeciesContainer species;

        string manifestFilename(argv[argc-1]);
        string dictFilename = manifestFilename + ".dict";

        ifstream ifs(manifestFilename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + manifestFilename);

        while (ifs) {
                string line;
                getline(ifs, line);
                if (line.empty())
                        continue;

                string speciesName, fastaFilename;
                istringstream iss(line);
                iss >> speciesName >> fastaFilename;

                if (speciesName.empty() || fastaFilename.empty()) {
                        throw runtime_error("File " + manifestFilename + " has incorrect format.\n"
                                            "Refer to the documention for more information.");
                }

                ifstream fn(fastaFilename);
                if (!fn)
                        throw runtime_error("Could not open fasta file: " + fastaFilename);

                species.addSeqFile(speciesName, fastaFilename);
        }

        ifs.close();

        // create the dictionary of species
        species.computeNucleotideComposition();
        species.write(dictFilename);

        cout << "Done!  Dictionary written to " << dictFilename << endl;
}
