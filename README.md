# blamm - BLAS Accelerated Motif Matching

BLAS Accelerated Motif Matching or blamm is a tool to find Position Weight Matrix (PWM) occurrences in DNA sequences. It runs on both CPU and GPU, requires very little RAM and can make using of multithreading. The software is written in C++.

## Getting Started

The performance of blamm depends almost exclusively on the performance of the BLAS library. BLAS (Basic Linear Algebra Subprograms) provide routines for basic vector and matrix operations. The blamm software relies heavily on matrix-matrix multiplications (`sgemm`).  Almost all CPU/GPU vendors provide highly optimized BLAS implementations.

* For Intel CPUs we highly recommend the use of MKL (Math Kernel Library).
* For AMD CPUs we highly recommend the use of ACML (AMD Core Math Library).
* The open source GotoBLAS or OpenBLAS implementations may also provide good performance.
* For GPUs, we support nVidia's CUBLAS.

We do *not* recommend the use of BLAS libraries as provided by many Linux distributions (often called `blas` and `blas-devel`). These reference implementations are often tens of times slower than optimized implementations.

### Prerequisites

* CMake version 2.80 or higher
* A C++ compiler with C++11 support (e.g. gcc version 4.8.1 or higher)
* BLAS library
* Optionally for the GPU version: nVidia CUDA SDK (we tested version 8.0)

### Installing 

Installing may be as easy as:

```bash
tar -xzvf blamm-X.Y.Z.tar.gz
cd blamm-X.Y.Z
mkdir build
cd build
cmake ..
make
```

The `blamm` binary file will then be present in the `blamm.X.Y.Z/build` directory.  Optinally, to install the binary system-wide:

```bash
make install
```

CMake will automatically detect the presence or absence of the C++ compiler, BLAS libraries and CUDA compiler. It will generate warnings/errors in case software is missing from your system.

### Running the example

Example input files are available in the `blam.X.Y.Z/example` folder.

```bash
# Generate the dictionary file for the input sequence(s)
blamm dict sequences.mf

# Generate the PWM histograms
blamm hist motifs.jaspar sequences.mf

# Find the PWM occurrences in the sequence(s)
blamm scan -rc -pt 0.0001 motifs.jaspar sequences.mf
```

Also available in that folder is the default `settings.cnf` file. This file contains certain advanced parameters to tweak performance. For most users, the default settings are fine.

## Manual

The `blamm` binary consists of three modules:

* `blamm dict` - generate a sequence dictionary and compute the background nucleotide composition
* `blamm hist` - generate PWM score histograms (optionally)
* `blamm scan` - find PWM occurrences in the sequence

For each module, help messages can be displayed using the `-h` or `--help` flag, e.g.:

```bash
blamm dict --help
```

By default, `blamm` uses as many threads as the operating system suggests (typically the number of CPU cores or twice that number in case the CPU has hyperthreading enabled). This may be changed using the `-t <val>` or `--numthreads <val>` option in the hist and scan modules.

### Dictionary module

The input sequence(s) must be prepared as one or more *fasta* files. Every sequence must be preceded by a sequence descriptor line that start with a `>` character. The sequence descriptors are used to report PWM matches. Hence, in order to avoid ambiguities, each sequence should have unique sequence descriptor. Note that sequence descriptors are read by `blamm` up until the first whitespace character. The remainder of the descriptor line ignored. For example, see this minimal [example](./example/sequences/fileA.fasta) fasta file.

To facilitate command line interaction with `blamm`, the input fasta files should be listed in a single manifest file. For example, an [example `sequences.mf` file could look like this:

```bash
group1		pathTo/fileA.fa
group2		pathTo/fileB.fa
group1		pathTo/fileC.fa
```

Where each fasta file is preceded by a user-defined *group identifier* (see further). One can now run the `blamm dict` module:

```bash
blam dict sequences.mf
```

`blamm` will read through all input fasta files and create a dictionary file with all sequence descriptors and the nucleotide (ACGT) compositions. The latter serve as background models b<sub>i</sub> when computing the Position Weight Matrices (PWMs). This procedure is light-weight and should complete in a few seconds. The results are written to the file `sequences.mf.dict` and will be used by the other modules.

Separate background nucleotide compositions are computed across all fasta files that share the same *group identifier* in the manifest file. In the example above, one background model is computed across all sequences contained in `fileA.fa` and `fileC.fa` as both files share the same group identifier. A second background model is computed for all sequences contained in `fileB.fa`. In case a single background model should be computed across all sequences in the three files, it suffices to provide them with the same group identifier in the manifest file.

Note that different sequences in the same fasta file will always use the same background model. In case the use of different background models is warranted, the sequences should be organised in different fasta files (each with a unique group identifier).

Note that the choice of the actual name of the group identifier(s) is arbitrary. The only restriction is that group identifiers should not contain whitespace characters.

### Hist module (optional)

The motif search patterns should be provided as Position Frequency Matrices (PFM) in JASPAR format.  All PFMs should be organized in a single file that ends with the `.jaspar` suffix. For example, see this minimal [example motif](./example/motifs.jaspar) file.

Note that JASPAR motifs are typically provided as one separate `.jaspar` file per motif. These can easily be merged using `cat *.jaspar > motifs.jaspar`.

The `blamm hist` module generates a PWM score histogram for each motif - background model combination. For example, if the motif input file contains three search patterns and two group identifiers (i.e., two background models) were specified in the manifest file, a total of six histograms will be computed. The histograms are important to convert p-values into PWM thresholds. Histograms can be computed as follows:

```bash
blamm hist motifs.jaspar sequences.mf
```
For each motif - background model combination, two files are generated: `hist_groupID_motifID.dat` and `hist_groupID_motifID.gnu`. The `.dat`file contains the histogram data while the `.gnu` file can be used by `gnuplot` to generated a postscript `.ps` file that can visualized by a viewer of your choise (e.g. `okular`):

```bash
gnuplot hist_groupID_motifID.gnu
okular hist_groupID_motif_ID.ps
```

In case no optional arguments are provided to the `blamm hist` module, theoretical histograms are computed using a dynamic programming algorithm. In that case, the background model nucleotide probabilities are used assuming a zero-order model.

Alternatively, one can generate *empirical* PWM score histograms by scanning a portion of the actual input sequences using the `-e` or `--empirical` flag. By default, the first 10 million nucleotides of the sequences in each group are scanned. This can be changed using the `-l <val>` or `--length <val>` flag.  For example, to scan the first 20 million nucleotides:

```bash
blamm hist -e -l 20000000 motifs.jaspar sequences.mf
```

Regardless of whether the *theorertical* or *empirical* histograms are computed, the user can specify the number of bins per histogram (250 by default) using the `-b <val>` or `--numbins <val>` flag. A higher number of bins improves the accuracy when converting p-values to PWM score thresholds.

Finally, the `-H <path>` or `--histdir <path>` option can be used to provide `blamm` with a path to an existing directory where the histogram output files can be written.

### Scan module

The scan module is used to scan the input sequences for PWM matches, i.e., substrings in the sequences for which the PWM score exceeds the PWM threshold. There are three options to specify the threshold:

* Using the `-at <T>` or `--absthreshold <T>` flag. This sets a global, absolute PWM score threshold T that is the same for all PWMs. This is generally not recommended as a *good* threshold differs between PWMs.
* Using the `-rt <f>` or `--relthreshold <f>` flag with f a fraction [0..1]. For each PWM individually, the score threshold is computed as (1-f) * minScore + f  * maxScore. As the shape of PWM score histograms may differ between PWMs, this choice is again not optimal.
* Using the `-pt <p>` or `--pthreshold <p>` flag with p the p-value. Using the histograms that were previously computed, the p-value will be converted to a PWM threshold for each PWM individually.

The scan module can be run as follows:

```bash
blamm scan -pt 0.00001 motifs.jaspar sequences.mf
```

The motif occurrences are written to `occurrences.txt`by default. This may be changed by providing the `-o <outputfilename>` or `--output <outputfilename>` options. Motif occurrences are written in the `.gtf` file format, example:

```bash
chr1    blamm   MA0282.1        22893   22901   13.5375 -       .       .
chr1    blamm   MA0282.1        29577   29585   12.8617 -       .       .
chr1    blamm   MA0282.1        57287   57295   11.453  -       .       .
chr1    blamm   MA0282.1        59932   59940   13.5375 +       .       .
...
```

The first column refers to sequence identifier (as specified in the input fasta file). The second column will always be blamm (the tool that generated the data). The third column is the motif identifier as specified in the motifs.jaspar input file. The fourth and fifth columns are respectively the begin and end positions motif occurrence in the sequence. The sixth column is the PWM score. Finally, the seventh column is the strand ('+' for forward; '-' for reverse). Column eight and nine are not used.

The PWM thresholds that were used are written to `PWMthresholds.txt`. This file is formatted as follows:

```bash
group1    MA0282.1        -37.9454        10.9633 13.5375
...
```

The first column refers to the group identifier (i.e., the background model). The second column contains the motif identifier. The third, fourth and firth column respectively contain the minimum, threshold and maximum PWM score.


In case the histograms were written to a different directory, the `-H <path>` or `--histdir <path>` option should again be used to point to their correct location. Note that the histograms are only needed in case a p-value threshold is specified.

By default, only the forward strand of the sequences is scanned. Using the `-rc` or `--revcompl` flag enables the search in both forward and reverse strands.

Finally, the `-c` or `--cuda` flag enables the use of GPU devices. By default, all available GPU devices will be used. By specifying the correct environment variables, the use of specific devices can be specified. We refer to nVidia documentation for this.

### Advanced settings

If the file is available, `blamm` reads the configuration file `settings.cnf` that is used to tweak certain performance parameters. If `settings.cnf` is not present, default values will be used. The file `settings.cnf` can be used to specify the following parameters:

```bash
MATRIX_S_W      250
MATRIX_S_H      1000
MATRIX_P_TILE_MIN_ZERO_AREA          4096
PSEUDOCOUNT     0.25
FLUSHOUTPUT     100000
```

PSEUDOCOUNT is the only parameters that may impact the reported motif occurrences by `blamm`. It corresponds to the value that is added to every element individually in the Position Frequency Matrix (PFM) before conversion to the Position Weight Matrix (PWM).

MATRIX_S_W and MATRIX_S_H specify the S-matrix dimensions. In the paper, these values are referred to as w and h respectively.  The MATRIX_P_TILE_MIN_ZERO_AREA specifies the minimum number of zeros that should be removed by splitting the P matrix in two tiles. FLUSHOUTPUT hints `blamm` how many motif occurrences should be kept in memory before spilling them to disk.


## Author

Jan Fostier, Ghent University - imec

## License

This project is licensed under the GNU General Public License Version 3 - see the [LICENSE](LICENSE) file for details.