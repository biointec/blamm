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

#include <string>
#include <iostream>

#include "dict.h"
#include "blstools.h"
#include "bls.h"
#include "pwmscan.h"
#include "ortho.h"
#include "hist.h"

using namespace std;

void printProgramVersion()
{
        cout << "BLAS Accelerated Motif Matching (blamm) -- version "
             << BLAMM_MAJOR_VERSION << "." << BLAMM_MINOR_VERSION
             << "." << BLAMM_PATCH_LEVEL << "\n";

        cout << "Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)\n";
        cout << "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE.\n" << endl;
}

void printUsage()
{
        cout << "Usage: blamm command [options]\n\n";

        cout << " command\n";
        cout << "  dict\t\t\tmake a dictionary for the input sequences\n";
        cout << "  hist\t\t\tgenerate PWM score histograms\n";
        cout << "  scan\t\t\tscan for pwm occurrences\n\n";
        //cout << "  bls\t\t\tcompute the branch length score\n";
        //cout << "  ortho\t\t\tconvert ortho groups\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << "  -v\t--version\tdisplay version\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void runDictModule(int argc, char **argv)
{
        try {
                Dictionary dict(argc, argv);
        } catch (runtime_error& e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

void runHistModule(int argc, char **argv)
{
        try {
                Histogram hist(argc, argv);
        } catch (runtime_error& e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

void runBLSModule(int argc, char **argv)
{
        try {
                BLS bls(argc, argv);
        } catch (runtime_error& e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

void runScanModule(int argc, char **argv)
{
        try {
                PWMScan scan(argc, argv);
        } catch (runtime_error& e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

void runOrthoModule(int argc, char **argv)
{
        try {
                Ortho ortho(argc, argv);
        } catch (runtime_error& e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

int main(int argc, char **argv)
{
        Command command = Command::none;

        // parse first parameter
        if (argc > 1) {
                string arg(argv[1]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-v") || (arg == "--version")) {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (arg == "dict") {
                        command = Command::dict;
                } else if (arg == "hist") {
                        command = Command::hist;
                } else if (arg == "scan") {
                        command = Command::scan;
                } /*else if (arg == "bls") {
                        command = Command::bls;
                } else if (arg == "ortho") {
                        command = Command::ortho;
                }*/
        }

        switch (command)
        {
                case Command::none:
                        printUsage();
                        exit(EXIT_FAILURE);
                        break;
                case Command::dict:
                        runDictModule(argc, argv);
                        break;
                case Command::hist:
                        runHistModule(argc, argv);
                        break;
                case Command::bls:
                        runBLSModule(argc, argv);
                        break;
                case Command::scan:
                        runScanModule(argc, argv);
                        break;
                case Command::ortho:
                        runOrthoModule(argc, argv);
                        break;
        }

        cout << "Exiting... bye!" << endl;

        return EXIT_SUCCESS;
}
