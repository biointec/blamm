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
#include <iomanip>
#include <thread>
#include <algorithm>
#include <sstream>
#include <future>
#include <cstring>
#include <functional>

#include "pwmscan.h"
#include "sequence.h"
#include "species.h"

#ifdef HAVE_CUDA
        #include <cuda_runtime.h>
#endif

using namespace std;

extern void kernel_wrapper(float *R, int m, int n, float *threshold, int* occIdx, float* occScore, int *nOcc);

void PWMScan::printUsage() const
{
        cout << "Usage: blamm scan [options] motifs.input sequences.input\n";
        cout << "Goal: find PWM matches in sequences\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help message\n";
        cout << "  -s\t--simple\tenable simple (slower) scan algorithm\n";
        cout << "  -c\t--cuda\t\tenable the CUDA scan algorithm\n";
        cout << "  -rc\t--revcompl\talso search the reverse strand for occurrences\n\n";

        cout << " [options arg]\n";
        cout << "  -at\t--absthreshold\tset the minimal absolute score for a motif occurrence\n";
        cout << "  -rt\t--relthreshold\tset the minimal relative score [0..1] for a motif occurrence (default = 0.95)\n";
        cout << "  -pt\t--pthreshold\tcompute the motif score threshold from p-value [0..1] (default = 1E-4)\n";
        cout << "  -t\t--numthreads\tset the number of parallel threads [default = #cores]\n\n";

        cout << " [file_options]\n";
        cout << "  -H\t--histdir\tdirectory where the histogram file(s) are stored [default = .]\n";
        cout << "  -o\t--output\tfilename for the motif occurrences [default = occurences.txt]\n\n";

        cout << " File \"motifs.input\" should contain the motifs in Jaspar format\n";
        cout << "  (see documentation for specification)\n\n";

        cout << " File \"sequences.input\" should contain a list of input fasta files in the following format:\n";
        cout << "   speciesID_1\tspecies1_sequences.fasta\n";
        cout << "   speciesID_2\tspecies2_sequences.fasta\n";
        cout << "   speciesID_3\tspecies2_sequences.fasta\n";
        cout << "   ...\n";
        cout << " where speciesID_x is a user-defined identifier per species\n";
        cout << " Refer to the documentation for more information\n\n";

        cout << " Example:\n";
        cout << "  blamm scan -o occurences.txt motifs.input sequences.input\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void PWMScan::writeOccToDisk(const std::vector<MotifOccurrence>& occurrences)
{
        const MotifContainer& mc = motifContainer;
        const SpeciesContainer& sc = speciesContainer;

        // produce a string outside the mutex
        ostringstream oss;
        for (auto o : occurrences) {
                oss << sc[o.getSpeciesID()].getSeqName(o.getSequenceID()) << "\t"
                    << "blamm\t"
                    << mc[o.getMotifID()].getName() << "\t"
                    << o.getSequencePos() << "\t"
                    << o.getSequencePos() + mc[o.getMotifID()].size() << "\t"
                    << o.getScore() << "\t"
                    << o.getStrand() << "\t.\t.\n";
        }

        // dump the string to disk under mutex protection
        lock_guard<mutex> lock(myMutex);
        totMatches += occurrences.size();
        os << oss.str();
}

void PWMScan::extractOccurrences(const Matrix& R, size_t offset,
                                 size_t speciesID, SeqMatrix& sm,
                                 vector<MotifOccurrence>& motifOcc)
{
        for (size_t j = 0; j < R.nCols(); j++) {

                // get the motif information
                size_t motifIdx = motifContainer.getMotifIDAtCol(j);
                const Motif& m = motifContainer[motifIdx];
                const float threshold = m.getThreshold();

                for (size_t i = 0; i < sm.getNumOccRow(); i++) {

                        float thisScore = R(i,j);
                        if (thisScore < threshold)
                                continue;

                        // at this point an occurrence is found
                        SeqPos seqPos = sm.getSeqPos(i, offset);
                        size_t remSeqLen = sm.getRemainingSeqLen(i, offset);

                        if (m.size() > remSeqLen)
                                continue;

                        char strand = m.isRevCompl() ? '-' : '+';
                        motifOcc.push_back(MotifOccurrence(motifIdx, speciesID, seqPos.getSeqIndex(),
                                                           seqPos.getSeqPos(), strand, thisScore));
                }
        }
}

void PWMScan::extractOccurrences2(int LDR, const map<int, int>& offset_v,
                                  int *occIdx, float *occScore, size_t speciesID,
                                  SeqMatrix& sm, vector<MotifOccurrence>& motifOcc)
{
        int idx = 0;
        for (const auto& it: offset_v) {
                int offset = it.first;
                int thisOcc = it.second;

                for (int c = 0; c < thisOcc; c++, idx++) {
                        float thisScore = occScore[idx];
                        int i = occIdx[idx] % LDR;        // row in R
                        int j = occIdx[idx] / LDR;        // col in R

                        size_t motifIdx = motifContainer.getMotifIDAtCol(j);
                        const Motif& m = motifContainer[motifIdx];

                        SeqPos seqPos = sm.getSeqPos(i, offset);
                        size_t remSeqLen = sm.getRemainingSeqLen(i, offset);

                        if (m.size() > remSeqLen)
                                continue;

                        char strand = m.isRevCompl() ? '-' : '+';
                        motifOcc.push_back(MotifOccurrence(motifIdx, speciesID, seqPos.getSeqIndex(),
                                                           seqPos.getSeqPos(), strand, thisScore));
                }
        }
}

void PWMScan::scanThreadNaive(size_t speciesID,
                              ProgressIndicator& progInd,
                              FastaBatch& seqBatch)
{
        vector<MotifOccurrence> occurrences;

        // we're reading the fasta files in blocks that overlap
        SeqBlock block;
        size_t overlap = motifContainer.getMaxMotifLen() - 1;
        size_t payload = settings.matrix_S_h * settings.matrix_S_w;

        while (seqBatch.getNextOverlappingBlock(block, payload, overlap)) {
                for (size_t i = 0; i < min<size_t>(block.size(), payload); i++) {
                        for (size_t j = 0; j < motifContainer.size(); j++) {
                                const Motif& m = motifContainer[j];
                                float thisScore = m.getScore(block.substr(i, m.size()));
                                if (thisScore < m.getThreshold())
                                        continue;

                                SeqPos seqPos = block.getSeqPos(i);
                                size_t remSeqLen = block.getRemainingSeqLen(i);

                                if (m.size() > remSeqLen)
                                        continue;

                                char strand = m.isRevCompl() ? '-' : '+';
                                occurrences.push_back(MotifOccurrence(j, speciesID, seqPos.getSeqIndex(),
                                                                      seqPos.getSeqPos(), strand, thisScore));
                        }

                        if (occurrences.size() > settings.flushOutput) {
                                writeOccToDisk(occurrences);
                                occurrences.clear();
                        }
                }

                writeOccToDisk(occurrences);
                occurrences.clear();

                progInd.printMessage(seqBatch.getProgressPerc());
        }
}

void PWMScan::scanPWMNaive(size_t speciesID, FastaBatch& seqBatch)
{
        // start histogram threads
        ProgressIndicator progInd;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadNaive, this,
                                          speciesID, ref(progInd), ref(seqBatch));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        progInd.finalize();
}

void PWMScan::scanThreadBLAS(size_t speciesID,
                             ProgressIndicator& progInd,
                             FastaBatch& seqBatch)
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

        vector<MotifOccurrence> occurrences;
        while (sm.getNextSeqMatrix(seqBatch)) {
                for (size_t offset = 0; offset < w; offset++) {
                        for (size_t i = 0; i < matrixTiles.size(); i++)
                                p.A_array[i] = sm.getData() + 4*offset*p.LDA[i];

                        Matrix::sgemm_batch(p);
                        extractOccurrences(R, offset, speciesID, sm, occurrences);

                        if (occurrences.size() > settings.flushOutput) {
                                writeOccToDisk(occurrences);
                                occurrences.clear();
                        }
                }

                // write the occurrences to disk
                writeOccToDisk(occurrences);
                occurrences.clear();

                progInd.printMessage(seqBatch.getProgressPerc());
        }
}

void PWMScan::scanPWMBLAS(size_t speciesID,
                          FastaBatch& seqBatch)
{
        // start scan threads
        ProgressIndicator progInd;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadBLAS, this,
                                          speciesID, ref(progInd), ref(seqBatch));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        progInd.finalize();
}

#ifdef HAVE_CUDA
void PWMScan::scanThreadCUBLAS(int devID, size_t speciesID,
                               ProgressIndicator& progInd, FastaBatch& seqBatch)
{
        float *d_P = 0, *d_S = 0, *d_R = 0;
        float *d_threshold = 0, *d_occScore = 0;
        int *d_occIdx = 0, *d_nOcc = 0;

        size_t overlap = motifContainer.getMaxMotifLen() - 1;
        size_t w = settings.matrix_S_w;
        size_t h = settings.matrix_S_h;

        // set CUDA device for this thread
        cudaSetDevice(devID);

        // create a CUBLAS handle
        cublasHandle_t handle;
        cublasCreate(&handle);

        // pattern matrix
        const Matrix& P = motifContainer.getMatrix();
        if (cudaMalloc((void **)&d_P, P.nRows() * P.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for P\n");
        cublasSetVector(P.nRows() * P.nCols(), sizeof(float), P.getData(), 1, d_P, 1);

        // sequence matrix
        SeqMatrix sm(h, w, overlap);
        if (cudaMalloc((void **)&d_S, 4*(w + overlap) * h * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for S\n");

        // result matrix
        Matrix R(h, P.nCols());
        if (cudaMalloc((void **)&d_R, R.nRows() * R.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for R\n");

        // set the thresholds
        float *threshold = new float[R.nCols()];
        for (size_t i = 0; i < R.nCols(); i++) {
                size_t motifID = motifContainer.getMotifIDAtCol(i);
                const Motif& m = motifContainer[motifID];
                threshold[i] = m.getThreshold();
        }

        if (cudaMalloc((void **)&d_threshold, R.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for threshold\n");
        cublasSetVector(R.nCols(), sizeof(float), threshold, 1, d_threshold, 1);
        delete [] threshold;

        // data structures for the occurrences
        vector<MotifOccurrence> occurrences;
        float *occScore = new float[2 * R.nRows() * R.nCols()];
        if (cudaMalloc((void **)&d_occScore, 2 * R.nRows() * R.nCols() * sizeof(float)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for occurrence scores\n");

        int *occIdx = new int[2 * R.nRows() * R.nCols()];
        if (cudaMalloc((void **)&d_occIdx, 2 * R.nRows() * R.nCols() * sizeof(int)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for occurrence indices\n");

        int nOcc = 0;
        if (cudaMalloc((void **)&d_nOcc, sizeof(int)) != cudaSuccess)
                throw runtime_error("Cannot allocate memory on CUDA device for number of occurrences\n");
        cublasSetVector(1, sizeof(int), &nOcc, 1, d_nOcc, 1);

        // set the BLAS parameters
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
                p.B_array[i] = d_P + matrixTiles[i].colStart*p.LDB[i];
                p.C_array[i] = d_R + matrixTiles[i].colStart*p.LDC[i];
        }

        map<int, int> offset_v;

        vector<future<void> > outputTask(numThreadsPerDevice);
        size_t currOutputTask = 0;

        while (sm.getNextSeqMatrix(seqBatch)) {
                // copy the sequence matrix to the device
                cublasSetVector(4*(w + overlap)*h, sizeof(float), sm.S.getData(), 1, d_S, 1);

                for (size_t offset = 0; offset < w; offset++) {
                        for (size_t i = 0; i < matrixTiles.size(); i++)
                                p.A_array[i] = d_S + 4*offset*p.LDA[i];

                        Matrix::sgemm_batch_cuda(handle, p);
                        kernel_wrapper(d_R, R.nRows(), R.nCols(), d_threshold,
                                       d_occIdx, d_occScore, d_nOcc);
                        int prevOcc = nOcc;
                        cublasGetVector(1, sizeof(int), d_nOcc, 1, &nOcc, 1);
                        offset_v[offset] = nOcc - prevOcc;

                        // don't retreive vector yet if insufficient results are on the GPU
                        if ( (nOcc <= int(R.nRows() * R.nCols())) && (offset < w-1) )
                                continue;

                        // get the results
                        cublasGetVector(nOcc, sizeof(float), d_occScore, 1, occScore, 1);
                        cublasGetVector(nOcc, sizeof(int), d_occIdx, 1, occIdx, 1);
                        extractOccurrences2(R.nRows(), offset_v, occIdx, occScore, speciesID, sm, occurrences);
                        offset_v.clear();
                        nOcc = 0;
                        cublasSetVector(1, sizeof(int), &nOcc, 1, d_nOcc, 1);
                }

                // write the output to disk
                if (outputTask[currOutputTask].valid())
                        outputTask[currOutputTask].get();
                outputTask[currOutputTask] = async(launch::async, &PWMScan::writeOccToDiskCopy, this, occurrences);
                currOutputTask = (currOutputTask + 1) % numThreadsPerDevice;

                occurrences.clear();

                progInd.printMessage(seqBatch.getProgressPerc());
        }

        // wait for all output to be written
        for (size_t i = 0; i < outputTask.size(); i++)
                if (outputTask[i].valid())
                        outputTask[i].get();

        delete [] occIdx;
        delete [] occScore;

        cudaFree(d_P);
        cudaFree(d_S);
        cudaFree(d_R);
        cudaFree(d_threshold);
        cudaFree(d_occScore);
        cudaFree(d_occIdx);
        cudaFree(d_nOcc);

        cublasDestroy(handle);
}

void PWMScan::scanPWMCUBLAS(size_t speciesID,
                            FastaBatch& seqBatch)
{
        // start one thread per GPU device
        ProgressIndicator progInd;
        vector<thread> workerThreads(numDevices);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&PWMScan::scanThreadCUBLAS, this, i,
                                          speciesID, ref(progInd), ref(seqBatch));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        progInd.finalize();
}
#endif

PWMScan::PWMScan(int argc, char ** argv) : simpleMode(false), cudaMode(false),
        totMatches(0), outputFilename("occurrences.txt"),
        absThSpecified(false), absThreshold(0.0), relThSpecified(false),
        relThreshold(0.95), pvalueSpecified(false), pvalue(0.0001),
        numThreads(thread::hardware_concurrency()), revCompl(false)
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
                } else if ((arg == "-rc") || (arg == "--revcompl")) {
                        revCompl = true;
                } else if ((arg == "-s") || (arg == "--simple")) {
                        simpleMode = true;
                } else if ((arg == "-c") || (arg == "--cuda")) {
                        cudaMode = true;
                } else if (((arg == "-at") || (arg == "--absthreshold")) && (i+1 < argc-2)) {
                        absThSpecified = true;
                        absThreshold = atof(argv[i+1]);
                        i++;
                } else if (((arg == "-rt") || (arg == "--relthreshold")) && (i+1 < argc-2)) {
                        relThSpecified = true;
                        relThreshold = atof(argv[i+1]);
                        if ((relThreshold < 0.0) || (relThreshold > 1.0))
                                throw runtime_error("The relative threshold should be in range [0..1].");
                        i++;
                } else if (((arg == "-pt") || (arg == "--pthreshold")) && (i+1 < argc-2)) {
                        pvalueSpecified = true;
                        pvalue = atof(argv[i+1]);
                        if ((pvalue < 0.0) || (pvalue > 1.0))
                                throw runtime_error("The p-value should be in range [0..1].");
                        i++;
                } else if (((arg == "-t") || (arg == "--numthreads")) && (i+1 < argc-2)) {
                        numThreads = atoi(argv[i+1]);
                        if (numThreads < 1)
                                throw runtime_error("Number of threads must be a non-zero positive number");
                        i++;
                } else if (((arg == "-H") || (arg == "--histdir")) && (i+1 < argc-2)) {
                        histdir = string(argv[i+1]);
                        if (histdir.back() != '/')
                                histdir.push_back('/');
                        i++;
                } else if (((arg == "-o") || (arg == "--output")) && (i+1 < argc-2)) {
                        outputFilename = string(argv[i+1]);
                        i++;
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        // If none of the thresholds is specified, default to relative threshold
        if (!(absThSpecified || relThSpecified || pvalueSpecified))
                relThSpecified = true;

        // Don't specify multiple thresholds
        if (absThSpecified && relThSpecified)
                throw runtime_error("Specify either the absolute or relative threshold, not both.");
        if (absThSpecified && pvalueSpecified)
                throw runtime_error("Specify either the absolute or p-value threshold, not both.");
        if (relThSpecified && pvalueSpecified)
                throw runtime_error("Specify either the relative or p-value threshold, not both.");

        cout << "Welcome to blamm -- PWM scan module" << endl;

        settings.printSettings();

        // A) load the manifest file
        string manifestFilename(argv[argc-1]);
        string dictFilename = manifestFilename + ".dict";
        speciesContainer.load(dictFilename);

        // B) load the motifs
        string motifFilename = string(argv[argc-2]);
        motifContainer.load(motifFilename, true);
        cout << "Loaded " << motifContainer.size() << " motifs from disk\n";
        cout << "Maximum motif size: " << motifContainer.getMaxMotifLen() << endl;

        if (revCompl) {
                motifContainer.addReverseComplements();
                cout << "Scanning both forward and reverse strand of the input sequence(s)" << endl;
        } else {
                cout << "Scanning only the forward strand of the input sequence(s)" << endl;
        }

        // compute the matrix tiles
        motifContainer.generateMatrixTiles(settings.matrix_P_tile_min_zero_area);

        // write some information about the threshold
        if (absThSpecified)
                cout << "Absolute motif score threshold set to: " << absThreshold << endl;
        else if (relThSpecified)
                cout << "Relative motif score threshold set to: " << relThreshold << endl;
        else if (pvalueSpecified)
                cout << "P-value motif score threshold set to: " << pvalue << endl;

        if (cudaMode) {
#ifdef HAVE_CUDA
                cudaGetDeviceCount(&numDevices);

                if (numDevices == 0) {
                        cerr << "CUDA error: no devices found. Aborting..." << endl;
                return;
                }

                cout << "Using " << numDevices << " GPU devices" << endl;

                numThreadsPerDevice = max<int>(numThreads / numDevices, 1);
                cout << "Using " << numThreadsPerDevice << " CPU threads per GPU device" << endl;
#else
                cerr << "ERROR: CUDA support not enabled\n";
                cerr << "Please recompile with CUDA support enabled" << endl;
#endif
        } else {
                cout << "Using " << numThreads << " thread(s)" << endl;
        }

        // C) scan the sequences
        ofstream ofsCutoff("PWMthresholds.txt");

        // start the output thread
        os.open(outputFilename.c_str());

        size_t speciesID = 0;
        for (auto species : speciesContainer) {
                cout << "Scanning species: " << species.getName();
                species.writeNuclProbabilities(settings.pseudocount);

                // compute PWMs and generate the pattern matrix
                motifContainer.generateMatrix(species.getNuclCounts(), settings.pseudocount);

                //motifContainer.writePossumFile("jaspar.possum");

                // set or compute the thresholds
                if (absThSpecified) {
                        for (auto& motif : motifContainer) {
                                motif.setThreshold(absThreshold);
                        }
                } else if (relThSpecified) {
                        for (auto& motif : motifContainer) {
                                float maxScore = motif.getMaxScore();
                                float minScore = motif.getMinScore();
                                float threshold = relThreshold * (maxScore - minScore) + minScore;
                                motif.setThreshold(threshold);
                        }
                } else if (pvalueSpecified) {
                        for (auto& motif : motifContainer) {
                                ScoreHistogram hist;
                                hist.loadHistogram(histdir, "hist_" + species.getName() + "_" +  motif.getBaseName());
                                motif.setThreshold(hist.getScoreCutoff(pvalue));
                        }
                }

                // write the thresholds to disk
                for (auto& motif : motifContainer) {
                        if (motif.size() > 15)
                                continue;
                        ofsCutoff << species.getName() << "\t" << motif.getName() << "\t"
                                  << motif.getMinScore() << "\t" << motif.getThreshold()
                                  << "\t" << motif.getMaxScore() << endl;
                }

                vector<string> filenames = species.getSequenceFilenames();
                FastaBatch seqBatch(filenames, species.getTotalSeqLength());

                // scan the sequences for PWM occurrences
                if (simpleMode) {
                        scanPWMNaive(speciesID++, seqBatch);
                } else if (cudaMode) {
#ifdef HAVE_CUDA
                        scanPWMCUBLAS(speciesID++, seqBatch);
#endif
                } else {
                        scanPWMBLAS(speciesID++, seqBatch);
                }
        }

        // close the output file
        os.close();

        ofsCutoff.close();

        cout << "\nWrote " << totMatches << " matches to " << outputFilename << ".\n";
}
