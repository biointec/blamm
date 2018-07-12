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
#include <sstream>
#include <stdexcept>

#include "bls.h"
#include "phylotree.h"

using namespace std;

void BLS::printUsage() const
{
        cout << "Usage: blstools bls [options] tree.newick leafs.input\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n\n";
        cout << "  -p\t--print\t\toutput tree to the screen\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\tfilename for output BLS values [default = stdout]\n\n";

        cout << " File \"tree.newick\" should contain the phylogenetic tree in Newick format, e.g.:\n";
        cout << "  \"(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);\"\n\n";

        cout << " File \"leafs.input\" should contain a list of space or tab separated\n";
        cout << "  leaf names for which the BLS should be computed (one line per entry)\n\n";

        cout << " Example:\n";
        cout << "  blstools bls -o output.bls tree.newick leafs.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

BLS::BLS(int argc, char** argv) : printTree(false)
{
        // check for sufficient arguments
        if (argc < 4) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // process optional arguments
        for (int i = 2; i < argc-2; i++) {
                string arg(argv[i]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-p") || (arg == "--print")) {
                        printTree = true;
                } else if (((arg == "-o") || (arg == "--output")) && (i+1 < argc-2)) {
                        outputFilename = string(argv[i+1]);
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        // process the input tree
        ifstream ifs(argv[argc-2]);
        if (!ifs)
                throw runtime_error("Could not open file: " + string(argv[argc-2]));
        string newickStr;
        getline(ifs, newickStr);
        ifs.close();

        PhylogeneticTree pt(newickStr);
        pt.normalizeBranchLength();

        if (printTree)
                cout << pt << endl;

        // open an output file if necessary
        ofstream ofs;
        if (!outputFilename.empty()) {
                ofs.open(outputFilename.c_str());
                cout.rdbuf(ofs.rdbuf());
        }

        // process the input leafs file
        ifs.open(argv[argc-1]);

        string line;
        while (getline(ifs, line)) {

                istringstream oss(line);
                set<string> leafs;

                string temp;
                while (oss >> temp)
                        leafs.insert(temp);

                cout << pt.getBLS(leafs) << endl;
        }

        ifs.close();

        if (!outputFilename.empty())
                ofs.close();
}
