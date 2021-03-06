# RibosePreferenceAnalysis
Preference analysis for ribose-seq data

This script is used for preference analysis of rNMP incorporation data generated from _ribose-seq_ protocol and Ribose-Map software. The neighbor of incorporated rNMP keeps certain preference in a specific species. This is because the surrounding dNMPs could affect the probability of rNMP misincorporation by replicative polymerases. This software is designed to reveal this preference.

## Citation
Balachander, S., Gombolay, A. L., Yang, T., Xu, P., Newnam, G., Keskin, H., El-Sayed, W., Bryksin, A. V., Tao, S., Bowen, N. E., Schinazi, R. F., Kim, B., Koh, K. D., Vannberg, F. O., & Storici, F. (2020). Ribonucleotide incorporation in yeast genomic DNA shows preference for cytosine and guanosine preceded by deoxyadenosine. _Nature communications, 11_(1), 2447. https://doi.org/10.1038/s41467-020-16152-5

## Dependency

Except Python3 standard libraries, the following packages are needed to run the scripts:

- Matplotlib

- Seaborn

- Scipy

## Usage

1. Use __dinuc_count.py__ to calculate background frequencies for the reference genome (__FASTA__ file).
   ```bash
   ./dinuc_count.py <ref genome> -o <background frequency>
   ```
   Dinucleotide background frequency is calculated by default. Available parameters are:
   1. __-d D__  Distance between dinucleotide pair, default = 1
   1. __-s__  Count only one strand
   1. __--mono__  Count mono nucleotide instead
   1. __--trinuc__  Count trinucleotide instead

1. Use __get_chrom.py__ to get background frequency of mitochondrial and nuclear DNA seperately
   ```bash
   ./get_chrom.py <background frequency> -s <chrM name> -o <chrM_frequency>
   ./get_chrom.py <background frequency> -s <chrM name> -v -o <nuc_frequency>
   ```
   Dinucleotide background frequency is calculated by default. Available parameters are:
   1. __-s S [S ...]__  Chromosome to be selected, default = chrM
   1. __-a__  Append to original file
   1. __--name NAME__  Name for the output line, default = input file name
   
1. Count different patterns of rNMP incorporation from __BED__ file generated by Ribose-Map using __dinuc_cnt_ribose.py__.
   ```bash
   ./dinuc_cnt_ribose.py <ref genome> <BED> -o <rNMP incorporation raw>
   ```
   By default only dinucleotide frequency is counted. Available parameters:
   1. __-f__  Use fourth column of bed file as frequency. Otherwise each row of BED file is considered as a rNMP
   1. __-m__  Also count mononucleotide frequency
   1. __-d__  Also count dinucleotide frequency
   1. __--dist DIST [DIST ...]__  Distance between rNMP and its dNMP neighbor
   1. __-t__  Also count trinucleotide frequency
   
1. Get the data of desired chromosome.
   ```bash
   ./get_chrom.py <infile1> <infile2> ... -s <chromosome1> <chromsome2> ... -o <file desired>
   ```
   Available parameters:
   1. __-v__  Select non-matching chromosomes 
   1. __-a__  Append to original file 
   1. __--name NAME__  Name for the output line, default = input file name

1. Normalization
   ```bash
   ./normalize.py <file desired> <chrM or nuclear frequency> -o <normalized file>
   ```
   Available parameters:
   1. __--group_len {0,4,16}__  Number of rows of which the sum is 1. If 0 is selected, the sum of all rows will be 1. default = 0.
   1. __--name NAME__  Name of chromosome in background frequency used for normalization, default = saccer
   
1. Rename and sort libraries if needed
   ```bash
   ./resort.py <normalized file> <order file> -o <sorted normalized file>
   ```
   Available parameters:
   1. __-d D__  Connector of library informations, default = '-'
   2. __-c C__  Column number of library name, default=1

2. Draw preference heatmaps.
   ```bash
   ./draw_ribose.py <normalized file> -o <output basename>
   ```
   A heatmap in __PNG__ format is generated by default. Available parameters:
   1. __-b B__  Select background file. If a file is selected, the background percentage is added to labels.
   1. __--nr__  Input file is NR type dinucleotide frequency. Of which the second position is rNMP.
   1. __--mono__  Input file is mononucleotide frequency.
   1. __--tri TRI__  Input file is trinucleotide frequency, with rNMP incorporated at TRI position
   1. __--background_chrom BACKGROUND_CHROM__  Chromosome name of background file, default = chrM
   1. __--legend_group {0,4,16}__  Number of lables of which the sum is 1. If 0 is selected, the sum of all labels will be 1. default = 0.
   1. __--no_annot__  Hide percentage annotation in each cell.
   1. __--cmax CMAX__  Maximum value in color scale. Any preferency beyond that will show as the maximum color.
   
## License

This software is under GNU GPL v3.0 license

## Contact

If you have any question, please contact me at [pxu64@gatech.edu](mailto:pxu64@gatech.edu).

