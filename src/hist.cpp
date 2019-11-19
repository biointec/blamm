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
#include <thread>
#include <algorithm>
#include <cmath>
#include <functional>

#include "hist.h"
#include "motif.h"
#include "sequence.h"
#include "matrix.h"
#include "species.h"

using namespace std;

void Histogram::printUsage() const
{
        cout << "Usage: blamm hist [options] motifs.input sequences.input\n";
        cout << "Goal: compute PWM score histograms\n\n";

        cout << " [options]\n";
        cout << "  -e\t--empirical\tcompute empirical histograms\n";
        cout << "  -h\t--help\t\tdisplay help message\n\n";

        cout << " [options arg]\n";
        cout << "  -b\t--numbins\tnumber of bins per histogram [default = 250]\n";
        cout << "  -t\t--numthreads\tset the number of parallel threads [default = #cores]\n";
        cout << "  -l\t--length\tlength of sequence data to analyze [default = 10000000]\n\n";

        cout << " [file_options]\n";
        cout << "  -H\t--histdir\toutput directory for the histogram file(s) [default = .]\n\n";

        cout << " File \"motifs.input\" should contain all motifs\n";
        cout << "  (see documentation for specification)\n\n";

        cout << " File \"sequences.input\" should contain a list of input"
             << " fasta files in the following format:\n";
        cout << "   speciesID_1\tspecies1_sequences.fasta\n";
        cout << "   speciesID_2\tspecies2_sequences.fasta\n";
        cout << "   speciesID_3\tspecies2_sequences.fasta\n";
        cout << "   ...\n";
        cout << " where speciesID_x is a user-defined identifier per species\n";
        cout << " Refer to the documentation for more information\n\n";

        cout << " Example:\n";
        cout << "  blamm hist motifs.input sequences.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void Histogram::extractObsScore(const Matrix& R, size_t offset,
                          const SeqMatrix& sm,
                          const MotifContainer& motifContainer,
                          vector<ScoreHistogram>& histContainer)
{
        for (size_t j = 0; j < R.nCols(); j++) {

                // get the motif information
                size_t motifIdx = motifContainer.getMotifIDAtCol(j);
                const Motif& m = motifContainer[motifIdx];

                for (size_t i = 0; i < sm.getNumOccRow(); i++) {
                        float thisScore = R(i,j);

                        // at this point an occurrence is found
                        size_t remSeqLen = sm.getRemainingSeqLen(i, offset);

                        if (m.size() > remSeqLen)
                                continue;

                        histContainer[motifIdx].addObservation(thisScore);
                }
        }
}

void Histogram::histThread(const MotifContainer& motifContainer,
                           FastaBatch& seqBatch,
                           vector<ScoreHistogram>& histContainer)
{
        size_t overlap = motifContainer.getMaxMotifLen() - 1;
        size_t w = settings.matrix_S_w;
        size_t h = settings.matrix_S_h;

        // pattern matrix
        const Matrix& P = motifContainer.getMatrix();

        // sequence matrix
        SeqMatrix sm(h, w, overlap);

        // result matrix
        Matrix R(h, P.nCols());

        // sgemm batch parameters
        const auto matrixTiles = motifContainer.getMatrixTiles();
        SgemmBatchParams p(matrixTiles.size());

        for (size_t i = 0; i < matrixTiles.size(); i++) {
                p.m[i] = h;
                p.k[i] = matrixTiles[i].rowEnd;
                p.n[i] = matrixTiles[i].colEnd-matrixTiles[i].colStart;
                p.LDA[i] = h;
                p.LDB[i] = P.nRows();
                p.LDC[i] = h;
                p.alpha[i] = 1.0f;
                p.beta[i] = 0.0f;
                p.B_array[i] = P.getData() + matrixTiles[i].colStart*p.LDB[i];
                p.C_array[i] = R.getData() + matrixTiles[i].colStart*p.LDC[i];
        }

        while (sm.getNextSeqMatrix(seqBatch)) {
                for (size_t offset = 0; offset < w; offset++) {
                        for (size_t i = 0; i < matrixTiles.size(); i++)
                                p.A_array[i] = sm.getData() + 4*offset*p.LDA[i];

                        Matrix::sgemm_batch(p);
                        extractObsScore(R, offset, sm, motifContainer, histContainer);
                }

                cout << "."; cout.flush();
        }
}

void Histogram::generateEmpiricalHist(const Species& species,
                                      const MotifContainer& motifContainer,
                                      vector<ScoreHistogram>& histContainer)
{
        vector<string> filenames = species.getSequenceFilenames();
        FastaBatch seqBatch(filenames, maxLength);

        // start histogram threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&Histogram::histThread, this,
                                          cref(motifContainer), ref(seqBatch),
                                          ref(histContainer));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << endl;
}

void Histogram::generateTheoreticalHist(const Species& species,
                                        const MotifContainer& motifContainer,
                                        vector<ScoreHistogram>& histContainer)
{
        for (size_t i = 0; i < motifContainer.size(); i++) {
                const Motif& m = motifContainer[i];
                map<float, float> spectrum;     // < score, PDF >
                array<float, 4> background = species.getNuclProbabilities(settings.pseudocount);
                m.computeTheoreticalSpectrum(numBins, background, spectrum);

                for (const auto& it : spectrum)
                        histContainer[i].setNumObservations(it.first, maxLength * it.second);
        }
}

Histogram::Histogram(int argc, char ** argv) : maxLength(10000000), numBins(250),
        numThreads(thread::hardware_concurrency()), empirical(false)
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
                } else if (((arg == "-l") || (arg == "--length")) && (i+1 < argc-2)) {
                        maxLength = atoll(argv[i+1]);
                        i++;
                } else if (((arg == "-b") || (arg == "--numbins")) && (i+1 < argc-2)) {
                        numBins = atoll(argv[i+1]);
                        if (numBins < 2)
                                numBins = 2;
                        i++;
                } else if ((arg == "-e") || (arg == "--empirical")) {
                        empirical = true;
                } else if (((arg == "-t") || (arg == "--numthreads")) && (i+1 < argc-2)) {
                        numThreads = atoi(argv[i+1]);
                        i++;
                } else if (((arg == "-H") || (arg == "--histdir")) && (i+1 < argc-2)) {
                        histdir = string(argv[i+1]);
                        if (histdir.back() != '/')
                                histdir.push_back('/');
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        cout << "Welcome to blamm -- histogram module" << endl;

        // A) load the manifest file
        string manifestFilename(argv[argc-1]);
        string dictFilename = manifestFilename + ".dict";

        SpeciesContainer speciesContainer;
        speciesContainer.load(dictFilename);

        // B) Load the motifs
        string motifFilename = string(argv[argc-2]);

        MotifContainer motifContainer;
        motifContainer.load(motifFilename, false);    // no permutations!
        cout << "Loaded " << motifContainer.size() << " motifs from disk";
        cout << "\nMaximum motif size: " << motifContainer.getMaxMotifLen() << endl;

        // compute the matrix tiles (only for empirical)
        if (empirical)
                motifContainer.generateMatrixTiles(settings.matrix_P_tile_min_zero_area);

        cout << "Using " << numThreads << " thread(s)" << endl;

        // C) compute the histograms
        for (auto species : speciesContainer) {
                cout << "Generating histograms for species: " << species.getName();
                species.writeNuclProbabilities(settings.pseudocount);

                // compute PWMs and generate the pattern matrix
                motifContainer.generateMatrix(species.getNuclCounts(), settings.pseudocount);

                // now generate score histograms for each motif
                vector<ScoreHistogram> histContainer;
                for (size_t i = 0; i < motifContainer.size(); i++) {
                        const Motif& m = motifContainer[i];
                        histContainer.push_back(ScoreHistogram(m.getMinScore(), m.getMaxScore(), numBins));
                }

                if (empirical)
                        generateEmpiricalHist(species, motifContainer, histContainer);
                else
                        generateTheoreticalHist(species, motifContainer, histContainer);

                // write all histograms to file
                for (size_t i = 0; i < histContainer.size(); i++)
                        histContainer[i].writeGNUPlotFile(histdir,
                                "hist_" + species.getName() + "_" +  motifContainer[i].getName(),
                                motifContainer[i].getName() + " (" + species.getName() + ")");
        }
}
