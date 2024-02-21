# MiSeq data analysis:

## Load the data files

We can import all of the zipped files output from the MiSeq together. These are already demultiplexed into different samples, with two files for each sample from the forward (R1) and reverse (R2) run. Can load all the sequences into QIIME2 together by specifying a folder conaining all the fastq.gz files.

Note that the data format from Earlham was a folder for each sample (with Fwd and Rev read in each), so needed to copy them all into a single folder first.

```
cp -t analysis PKG\ -\ ICL.TS.ENQ-5291.A.01_\ 376\ x\ v3-v4\ amplicons\ R0971\ S0001-S0376/230221_M01013_0020_000000000-KWM37/*/*.fastq.gz
```

Also due to the file format Earlham sent, need to import into qiime as a manifest, which requires a manifest file (sample-manifest.tsv) which looks like this:

```
sample-id     forward-absolute-filepath       reverse-absolute-filepath
sample-1      $PWD/some/filepath/sample0_R1.fastq.gz  $PWD/some/filepath/sample1_R2.fastq.gz
sample-2      $PWD/some/filepath/sample2_R1.fastq.gz  $PWD/some/filepath/sample2_R2.fastq.gz
sample-3      $PWD/some/filepath/sample3_R1.fastq.gz  $PWD/some/filepath/sample3_R2.fastq.gz
sample-4      $PWD/some/filepath/sample4_R1.fastq.gz  $PWD/some/filepath/sample4_R2.fastq.gz
```

Load in the data:

```
conda  activate qiime2-2022.2

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path sample-manifest.tsv --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
```

## Sequence quality control

Begin by performing quality control on the demultiplexed sequences, do this with DADA2.

This is a summary to visualise:

```
qiime demux summarize --i-data paired-end-demux.qza --o-visualization demux-summary.qzv
```

Can open these visualisations at https://view.qiime2.org/

The quality looks great for the forward read, drops off a bit towards the end of the reverse reads, but should still be adequate for joining reads I think.

Run the DADA2 plug-in. It takes paired-end-demux.qza as input and requires two parameters, trunc-len-f and trunc-len-r. The idea is to optimize merging of the forward and reverse reads by removing as much of the lower quality portions of the reads as possible and still leave enough overlap. Doing this by inspection of the quality plots is subjective, can try to use Zymo Research’s program FIGARO to find the parameters for us. See tutorial on FIGARO for how to install and run Figaro: http://john-quensen.com/tutorials/figaro/

We'll truncate the reverse reads to a shorter length than the forward reads due to their diminished quality. Also need to chop a bit off the start of the sequences because the primers are still in there (17 bases for the fwd primer and 21 bases for the reverse). 

```
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza --p-trim-left-f 17 --p-trunc-len-f 280 --p-trim-left-r 21 --p-trunc-len-r 225 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --verbose
```

Generate qiime artifacts to view summary statistics:

```
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
```
Now turn it into a featuredata table that we can use later on. Note we need to have created a metadata file for this (metadata tutorial here: https://docs.qiime2.org/2022.8/tutorials/metadata/).

The `feature-table summarize` command will give information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics. The `feature-table tabulate-seqs` command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.

```
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
```

Rarefaction to retain 10,000 sequences per sample seems like a decent option - we lose a couple of samples, but retain the vast majority.

# Downstream analysis

Now we want to:

* Call taxonomy
* Build a phylogeny
* Calculate diversity metrics
* Probably some other stuff...

I'm going to do all of this using the Python API plugin for qiime.

```
# import some stuff, including qiime2 commands
from qiime2 import Metadata
from urllib import request
import qiime2.plugins.metadata.actions as metadata_actions
from qiime2.plugins import feature_table
from qiime2 import Artifact
import qiime2.plugins.feature_classifier.actions as feature_classifier_actions
import qiime2.plugins.taxa.actions as taxa_actions
import qiime2.plugins.feature_table.actions as feature_table_actions
import qiime2.plugins.diversity.actions as diversity_actions
```
## Taxonomy assignment

```
## get the taxonomic classifier
url = 'https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/data/030-tutorial-downstream/020-taxonomy/gg-13-8-99-nb-classifier.qza'
fn = '../gg-13-8-99-nb-classifier.qza'
request.urlretrieve(url, fn)
gg_13_8_99_nb_classifier = Artifact.load(fn)

# load the data
filtered_sequences = Artifact.load("rep-seqs.qza")
filtered_table = Artifact.load("table.qza")
sample_metadata_md = Metadata.load("../metadata.tsv")

# use the classifier to assign taxonomic information to the sequences
taxonomy, = feature_classifier_actions.classify_sklearn(
    classifier=gg_13_8_99_nb_classifier,
    reads=filtered_sequences,
)

# and generate human-readable summary
taxonomy_as_md_md = taxonomy.view(Metadata)
taxonomy_viz, = metadata_actions.tabulate(
    input=taxonomy_as_md_md,
)

# write to file
taxonomy.save('taxonomy-silva.qza')

# Visualise
taxa_bar_plots_1_viz = taxa_actions.barplot(
    table=filtered_table,
    taxonomy=taxonomy,
    metadata=sample_metadata_md,
)

taxa_bar_plots_1_viz.visualization.save('taxa_bar_plots.qzv')
```

This is not perfect and needs to be re-done later. We've sequenced the v3-v4 16S region, but this taxonomic classifier is only for the v4 region. I tried a classifier using full 16S sequences from SILVA, but alas, my computer ran out of memory trying to process it. Will need to try to get it working on the HPC. For now, consider the taxonomy probably fine at high taxonomic levels, but probably a bit incorrect at species level.

## Phylogenetic tree construction

Next we’ll build a phylogenetic tree from our sequences using the q2-phylogeny plugin’s align_to_tree_mafft_fasttree action. This action is a pipeline in QIIME 2, which means that it strings together multiple simpler operations that are often performed together to reduce the number of steps that users have to take. Pipelines are used in the same way as other actions.

This particular pipeline performs four distinct steps steps:

* perform a multiple sequence alignment using mafft
* filter highly variable positions from the alignment (these positions tend to introduce noise into the phylogenetic tree)
* build an unrooted phylogenetic tree
* add a root to the unrooted tree

The final unrooted phylogenetic tree will be used for analyses that we perform next - specifically for computing phylogenetically aware diversity metrics. Output artifacts will be available for each of these steps.

```
import qiime2.plugins.phylogeny.actions as phylogeny_actions

action_results = phylogeny_actions.align_to_tree_mafft_fasttree(
    sequences=filtered_sequences,
)
aligned_rep_seqs = action_results.alignment
masked_aligned_rep_seqs = action_results.masked_alignment
unrooted_tree = action_results.tree
rooted_tree = action_results.rooted_tree

# write out the artifacts for later
aligned_rep_seqs.save("alignment.qza")
masked_aligned_rep_seqs.save("masked_alignment.qza")
unrooted_tree.save("tree.qza")
rooted_tree.save("rooted_tree.qza")
```

## Sampling depth/rarefaction

As we begin performing more analyses of the samples in our feature table, an important parameter that needs to be defined is the even sampling (i.e. rarefaction) depth that diversity metrics need to be computed at. Because most diversity metrics are sensitive to different sampling depths across different samples, it is common to randomly subsample the counts from each sample to a specific value. For example, if you define your sampling depth as 500 sequences per sample, the counts in each sample will be subsampled without replacement so that each sample in the resulting table has a total count of 500. If the total count for any sample(s) are smaller than this value, those samples will be dropped from the downstream analyses. Choosing this value is tricky. Its recommend to choose by reviewing the information presented in the feature table summary file (generated this previously). Choose a value that is as high as possible (to retain more sequences per sample) while excluding as few samples as possible.

I'll start with a sampling depth of 10000, to retain most of our samples.

### Alpha rarefaction plots

After choosing an even sampling depth, it’s also helpful to see if your diversity metrics appear to have stabilized at that depth of coverage. You can do this for alpha diversity using an alpha rarefaction plot.

```
shannon_rarefaction_plot_viz, = diversity_actions.alpha_rarefaction(
    table=filtered_table,
    metrics={'shannon'},
    metadata=sample_metadata_md,
    max_depth=10000,
)

shannon_rarefaction_plot_viz.save("shannon_rarefaction_plot.qzv")
```

## Core phylogenetic diversity metrics

core-metrics-phylogenetic requires a feature table, a rooted phylogenetic tree, and the sample metadata as input. It additionally requires the sampling depth that this analysis will be performed at (see above).

```
action_results = diversity_actions.core_metrics_phylogenetic(
    phylogeny=rooted_tree,
    table=filtered_table,
    sampling_depth=10000,
    metadata=sample_metadata_md,
)
rarefied_table = action_results.rarefied_table
faith_pd_vector = action_results.faith_pd_vector
observed_features_vector = action_results.observed_features_vector
shannon_vector = action_results.shannon_vector
evenness_vector = action_results.evenness_vector
unweighted_unifrac_distance_matrix = action_results.unweighted_unifrac_distance_matrix
weighted_unifrac_distance_matrix = action_results.weighted_unifrac_distance_matrix
jaccard_distance_matrix = action_results.jaccard_distance_matrix
bray_curtis_distance_matrix = action_results.bray_curtis_distance_matrix
unweighted_unifrac_pcoa_results = action_results.unweighted_unifrac_pcoa_results
weighted_unifrac_pcoa_results = action_results.weighted_unifrac_pcoa_results
jaccard_pcoa_results = action_results.jaccard_pcoa_results
bray_curtis_pcoa_results = action_results.bray_curtis_pcoa_results
unweighted_unifrac_emperor_viz = action_results.unweighted_unifrac_emperor
weighted_unifrac_emperor_viz = action_results.weighted_unifrac_emperor
jaccard_emperor_viz = action_results.jaccard_emperor
bray_curtis_emperor_viz = action_results.bray_curtis_emperor

# save all the new features
rarefied_table.save("rarefied/rarefied_table.qza")
faith_pd_vector.save("rarefied/faith_pd.qza")
observed_features_vector.save("rarefied/observed_features.qza")
shannon_vector.save("rarefied/shannon.qza")
evenness_vector.save("rarefied/evenness.qza")
unweighted_unifrac_distance_matrix.save("rarefied/unweighted_unifrac.qza")
weighted_unifrac_distance_matrix.save("rarefied/weighted_unifrac.qza")
jaccard_distance_matrix.save("rarefied/jaccard.qza")
bray_curtis_distance_matrix.save("rarefied/bray_curtis.qza")
unweighted_unifrac_pcoa_results.save("rarefied/unweighted_unifrac_pcoa.qza")
weighted_unifrac_pcoa_results.save("rarefied/weighted_unifrac_pcoa.qza")
jaccard_pcoa_results.save("rarefied/jaccard_pcoa.qza")
bray_curtis_pcoa_results.save("rarefied/bray_curtis_pcoa.qza")

# finally filter the sequences to this rarefied table too and save that
rarefied_sequences, = feature_table_actions.filter_seqs(
    data=filtered_sequences,
    table=rarefied_table,
)

rarefied_sequences.save("rarefied/rarefied_sequences.qza")
```


# Post-processing

Back into command-line qiime to output these as artifacts for use in picrust:

```
qiime tools export --input-path rarefied_table.qza --output-path rarefied_table
qiime tools export --input-path rarefied_sequences.qza --output-path rarefied_sequences
```

# PiCRUST Metagenome analysis

Now run those through picrust to try to assign function to sequences. This is in `data/Earlham/picrust`

First place the reads into a reference tree.
```
conda activate picrust2

place_seqs.py -s ../rarefied_sequences/dna-sequences.fasta -o out.tre -p 4 --intermediate intermediate/place_seqs -t sepp
```
Some sequences were discarded because they didn't align well to the reference sequences:

```
Warning - 266 input sequences aligned poorly to reference sequences (--min_align option specified a minimum proportion of 0.8 aligning to reference sequences). These input sequences will not be placed and will be excluded from downstream steps.
```

Then we do hidden-state prediction of gene families:
```
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n -e 0
hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1 -e 0
```
The output files of these commands are `marker_predicted_and_nsti.tsv.gz` and `EC_predicted.tsv.gz`.


`marker_predicted_and_nsti.tsv.gz`: The first column is the name of ASV, followed by the predicted number of 16S copies per ASVs, followed finally by the NSTI (nearest-sequenced taxon index) value per ASV. ASVs with a NSTI score above 2 are usually noise. It can be useful to take a look at the distribution of NSTI values for your ASVs to determine how well-characterized the community is overall and whether there are any outliers.

`EC_predicted.tsv.gz`: In this table the predicted copy number of all Enzyme Classification (EC) numbers is shown for each ASV. The NSTI values per ASV are not in this table since we did not specify the -n option. EC numbers are a type of gene family defined based on the chemical reactions they catalyze. For instance, EC:1.1.1.1 corresponds to alcohol dehydrogenase.

Now generate metagenome predictions. In the last step we predicted the copy numbers of gene families for each ASV. This output alone can be useful; however, often users are more interested in the predicted gene families weighted by the relative abundance of ASVs in their community. In other words, they are interested in inferring the metagenomes of the communities. This output can be produced by plugging in the BIOM table of ASV abundances per samples.
```
metagenome_pipeline.py -i ../rarefied_table/feature-table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out
```
  - `EC_metagenome_out/pred_metagenome_unstrat.tsv.gz` - overall EC number abundances per sample.
  - `EC_metagenome_out/pred_metagenome_contrib.tsv.gz` - A stratified table in "contributional" format breaking down how the ASVs contribute to gene family abundances in each sample.
  - `EC_metagenome_out/seqtab_norm.tsv.gz` - the ASV abundance table normalized by predicted 16S copy number.
  - `EC_metagenome_out/weighted_nsti.tsv.gz` - the mean NSTI value per sample (when taking into account the relative abundance of the ASVs). This file can be useful for identifying outlier samples in your dataset. In PICRUSt1 weighted NSTI values < 0.06 and > 0.15 were suggested as good and high, respectively. The cut-offs can be useful for getting a ball-park of how your samples compare to other datasets, but a weighted NSTI score > 0.15 does not necessarily mean that the predictions are meaningless.

The last major step of the PICRUSt2 pipeline is to infer pathway-level abundances. By default this script infers MetaCyc pathway abundances based on EC number abundances, although different gene families and pathways can also be optionally specified. This script performs a number of steps by default, which are meant to be as similar as possible to the approach implemented in HUMAnN2
```
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz -o pathways_out -p 1
```
Finally, it can be useful to have a description of each functional id in the output abundance tables. The below commands will add these descriptions as new column in gene family and pathway abundance tables.

```
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz
```

Great, now we can analyze this data in R.
