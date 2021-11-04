# Microbial populations are shaped by dispersal and recombination in a low biomass subseafloor habitat

Scripts used to process data and generate figures for the 2022 version of manuscript. For this repo, scripts have been modified to work with a limited number of examples. We apologize that the scripts are not more polished and are in slightly different styles. The code might not be useful for running your samples without modification, but there should be suitable information to customize scripts that process your data in a similar manner. Both Drs. Anderson and Tully are biologists before they are coders! :) 

##RPKM-plots

Dependencies: `numpy`, `pandas`, `seaborn`, `matplotlib`

`make_bubble_plots_binsofinterest.py` reads in the summary file `SNV-SAAV-bloomers-persisters-RPKMmin10-NEW-RPKM.txt` which contains multiple column headers but produces a plot of polymorphisms over time with data point size related to RPKM.

##allele-frequency-plots

Dependencies: `numpy`, `pandas`, `seaborn`, `matplotlib`

`make_allelefreq_clustermaps_allSNVs.py` reads in a anvi'o SNV table - in this instance the SNV table has been filtered to include positive entropy only. For the North Pond dataset, the input command explicitly targets the samples of interest. For example:

```
python make_allelefreq_clustermaps_allSNVs.py NORP83_SNV_allele_tracking_allSNVs.txt s82A_T1_SORTED_FILTERED,s82A_T3_SORTED_FILTERED,s82A_T4_SORTED_FILTERED,s82A_T6_SORTED_FILTERED
```

##gene-frequency-calc

Dependencies: `numpy`, `pandas`, `seaborn`, `matplotlib`

The `*_gene_cov_detection-GENE-COVERAGES.txt` is created as an output from anvi'o.

The script `get_gene_freq_from_anvio_output.py` converts file into the format displayed in `*_gene_cov_detection-GENE-COVERAGES-genefreq-script.txt`

```
python get_gene_freq_from_anvio_output.py *_gene_cov_detection-GENE-COVERAGES.txt
```

The script `systematic_genecov_ID.py` performs the math to identify genes that change coverage values >1x in the samples of interest. Sample IDs match the headerline. The output is a list stored as `*_gene_cov_detection-GENE-COVERAGES-genefreq-script.txt`

```
python systematic_genecov_ID.py *_gene_cov_detection-GENE-COVERAGES-genefreq-script.txt s82A_T0,s82A_T1,s82A_T2,s82A_T7,s82A_T8
```

The script `make_genecov_clustermap_specific_genes.py` creates a clustered heatmap of the genes of interest.

```
python make_genecov_clustermap_specific_genes.py *_gene_cov_detection-GENE-COVERAGES-genefreq-script.txt s82A_T0,s82A_T1,s82A_T2,s82A_T7,s82A_T8
```

##DESMAN-calc

Dependencies: `DESMAN`, `python==2.X`, `pandas`

DESMAN requires python2.X - i.e., scripts have old print statement formats.

Step 1: `get_variants_SCGs.py` - requires the list of Alenberg single copy COGs - `alneberg_single_copy_COGs.txt` - and then both the anvi'o generate COG table - `NORP83_COG_function.txt` - and the anvi'o generate SNV output, in this case filtered for positive entropy - `NORP83.snv.pos_entropy.txt`

```
python get_variants_SCGs.py NORP83_COG_function.txt NORP83.snv.pos_entropy.txt
```

Step 2: `convert_anvio_to_DESMAN.py` - this reformats the output created in Step 1 to be accessible by DESMAN

```
python convert_anvio_to_DESMAN.p NORP83.snv.pos_entropy_SCGonly.txt
```

Step 3: Reformat `NORP83.snv.pos_entropy_SCGonly_DESMANfriendly.txt` into appropriate DESMAN input files - with the format `*_outsel_var.csv` and `*_outtran_df.csv`

Step 4: In a directory for the MAG of interest - in this instanct NORP83 - run `run_DESMAN_triplicate.py`. The correctly formatted DESMAN file must be available in the directory.

```
python run_DESMAN_triplicate.py NORP83
```

Step 5: In a directory for the MAG of interest - `compile_posterior_dist_DESMAN.py` - compiles the posterior probabilities for the multiple iterations for further interpretation.

##inStrain-calc
Dependencies: `inStrain`

Requires access to filtered BAM files: `SORTED-FILTERED.tar.gz`. They can be accessed through the figshare repository: ADDRESS

## inStrain-related-calc
Dependencies: `pandas`, `numpy`, `tqdm`, `scikit-allel==1.3.3`, `biopython`, `scipy`, `seaborn`, `jupyter`

Both `calculate-four-gamete-test.py` and `calculate-Fst-from-inStrain-SNV-table.py` expect to access a directory called `/inStrain-Results` with multiple outputs representing different samples.

`Determine-Elevated-Fst-Window.ipynb` assumes that `*Fst.tsv` are created by `calculate-Fst-from-inStrain-SNV-table.py` prior to running.

##mcorr-calc
Dependencies: `mcorr`


Requires access to filtered BAM files: `SORTED-FILTERED.tar.gz`. They can be accessed through the figshare repository: ADDRESS

##gretel-calc
Dependencies: `gretel.yml` file contains expected dependencies

Gretel is a single threaded application that can be VERY time consuming for large datasets, this script uses `concurrent` to multithread the command, speeding up computation time

Requires access to filtered BAM files: `SORTED-FILTERED.tar.gz`. They can be accessed through the figshare repository: ADDRESS

##dNdS-calc
Dependencies: `biopython`, `concurrent`, `cd-hit`, `pal2nal.pl`, `muscle`, `codeml`, `plotly`

`batch-prep-dNdS.py` consumes the Gretel output directories and performs CD-HIT to determine the number of identical haplotypes and isoforms created. Then runs PAL2NAL and creates a CODEML control file.

A loop can be created to run through each gene using BASH:

```
ls -d dNdS-*/*/ > genedirectories.txt

while read p; do codeml "$p"codeml.ctl; mv rub 4fold.nuc 2ML.* rst* 2NG.* "$p"; done < genedirectories.txt
```

`make-dNdS-plots.py` reads in the `rst` file created by CODEML/PAML and filters the dN/dS results



