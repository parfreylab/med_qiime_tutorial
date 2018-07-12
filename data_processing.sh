#########################################
## MED & QIIME tutorial using lab data ##
#########################################

#### OVERVIEW ########

######
# Note: this tutorial begins after data has been filtered and validated
# For an overview of this process, please see steps 1-3 of this lab protocol: https://docs.google.com/document/d/1L83TgSxGLkJ3yyNPexzBB1qtrDlyrj1JwWiDEDjXusc/edit?usp=sharing
######
# This dummy dataset consists of 30 16S barcoded microbial samples from Westbeach, and 5000 reads were selected at random for each sample with seqtk to help make the MED runtime shorter
# We will be examining the differences between the communities

# In order to run many scripts, please invoke:
macqiime

# Demultiplex data: set quality threshold at 20
cd data
mkdir Westbeach
multiple_split_libraries_fastq.py -i raw_fastq -o Westbeach --demultiplexing_method sampleid_by_file

## PREPARING DATA FOR MED ##
# Trim data to 250bp and override default CCAATTGG
fastx_trimmer -l 250 -i seqs.westbeach.fna | fastx_clipper -v -a NNNNNNNNNNNNN -l 250 -o Westbeach/seqs.trimmed_clipped.fna > fastx_trim_clip.log
# To check number of reads after trimming:
wc -l seqs.trimmed_clipped.fna

# Modify headers by removing extra infofmation & format headers using trimheaders.rm_lead_zeroes.py in preparation for MED
# trimheaders.rm_lead_zeroes.py can be found in lab_scripts/trimheaders.rm_lead_zeroes.py on Github
# Feel free to read through the script to find out what the script does
cd ../..
python lab_scripts/trimheaders.rm_lead_zeroes.py -f data/Westbeach/seqs.trimmed_clipped.fna -o data/Westbeach/

###############
## USING MED ##
###############

# The first step after preparing our data is to run MED. 
# MED is invoked like this:
# Based on our number of reads, we picked -M 25 (most read sets require an -M value of at least 250):
decompose data/Westbeach/seqs.trimmed_clipped.MED.fna -M 25 -o MED/decompose-M_25
# To check number of OTU:
# Go into the MED/HTML-OUPUT/index.html output at each -M and check number of sequences represented. 
# Final OTU count: 310

# MED sorts sequences into groups (nodes) then calculates the entropy of the nodes based on sequence similarity. 
# each node is divided in two to reduce entropy, and the process is repeated until a minimum threshold of entropy is reached. 
# parameters that can be adjusted include:
# *minimum node size (how many sequences in entire node)
# *minimum node core size (how many sequences belong to node center)
# *maximum variation allowed per node (effects entropy threshold)

# when choosing a non-default M, look closely at your data, at the number of reads per each sample after you do your sequence trimming and filtering.
# generally, it's a good idea to make sure that your -M value won't eliminate taxa that are present in samples that have lower than average numbers of reads.

# example: I want to capture taxa that are only represent 1% of my total community population, and may only be present in a handful of samples. My data are variable in their read counts, some have over 50,000 reads but many only have 10,000 or less.
# I look at my summary of reads per sample, and choose an M value that corresponds to ~1% of the reads for one of my samples that has few reads but still has enough data for me to keep.

##############################
## PREPARING DATA FOR QIIME ##
##############################

# Steps for preparation:
# Transpose sample/OTU matrix
# MED adds a unique identifier for each node, but leaves leading zeroes in the IDs, so they are all the same character length (ex. 001, 058), we remove leading zeros before using QIIME
# MED also includes leading zeros in it's node representative fasta file, and the number of reads that belong to each node. QIIME doesn't like either of these things, so we remove them
# This can all be done by using the lab-coded script "trimheaders.rm_lead_zeroes.py"
# Note: this is the same script as used above
python lab_scripts/trimheaders.rm_lead_zeroes.py -m MED/decompose-M_25/MATRIX-COUNT.txt -n MED/decompose-M_25/NODE-REPRESENTATIVES.fasta -o data/Westbeach 

######################
## CHIMERA CHECKING ##
######################
## Method 1:
# 16s SILVA data can be found on the lab computers. Below are examples of a path to the files needed (example/path/), but please change as needed
# Check for any chimeric sequences using qiime and usearch61 for smaller datasets (in this example, assume data is 18s):
identify_chimeric_seqs.py -i data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.fasta -m usearch61 -o usearch_checked_chimeras/ -r example/path/../../../taxonomy_databases/SILVA_128/18s/rep_set_unaligned/99_SILVA_128_euks_rep_set.fasta

# Check number of sequences marked as chimeric. Usually no more than 10%
less usearch_checked_chimeras/identify_chimeric_seqs.log 

# Filter out the chimeric sequences 
filter_fasta.py -f data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.fasta -o data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras.fasta -s usearch_checked_chimeras/non_chimeras.txt

## Method 2:
# identify_chimeric_seqs.py doesn't work for large datasets (high memory footprint), so will need to use vsearch instead if you run out of memory
vsearch --uchime_ref data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.fasta --db example/path/../../../taxonomy_databases/SILVA_128/16s/rep_set_unaligned/99_SILVA_128_prks_rep_set.fasta --chimeras data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.chimeras.fasta --nonchimeras data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras.fasta
# Check number of sequences marked as chimeric. Usually no more than 10%
# In this case, there were 1.9% chimeras, so we can continue 

# Get list for filtering OTU table
grep ">" data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras.fasta | sed 's/>//g' > data/Westbeach/non_chimeras.txt 

# There's a log file from which you can check the output

########################
## ASSIGNING TAXONOMY ##
########################
# Make sure to use correct database (in this case, 16s)
# Assign taxonomy based on 16s taxonomy map and unaligned 16s sequences 
assign_taxonomy.py -i data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras.fasta -t example/path/../../../taxonomy_databases/SILVA_128/16s/taxonomy_maps/99_SILVA_128_taxa_map_7_levels.txt -r example/path/../../../taxonomy_databases/SILVA_128/16s/rep_set_unaligned/99_SILVA_128_prks_rep_set.fasta -m uclust -o ./assign_taxonomy

####################
# MAKING OTU TABLE #
####################
mkdir otu_table
biom convert -i data/Westbeach/MATRIX-COUNT.transposed.txt -o otu_table/OTU_Table.biom --to-json --table-type="OTU table" 

# Remove chimeric OTUs from the OTU table prior to diversity analyses.
filter_otus_from_otu_table.py -i otu_table/OTU_Table.biom -o otu_table/OTU_Table.no_chimeras.biom -e data/Westbeach/non_chimeras.txt --negate_ids_to_exclude

# Now check to see if any samples need to be filtered due to low reads, etc.:
biom summarize-table -i otu_table/OTU_Table.biom -o otu_table/OTU_Table_summary.txt
# Open above file. Did not observe any low reads (below 1000), so we will use all the data and have no need to filter
# Filtering based on hard-read threshold is up to the user

# Add taxonomy to OTU table
biom add-metadata -i otu_table/OTU_Table.biom -o otu_table/OTU_Table.wtaxa.biom --observation-metadata-fp assign_taxonomy/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy

############################
# MAKING PHYLOGENETIC TREE #
############################
# Using NODE-REPRESENTATIVES.DOWNSTREAM.fasta to assign taxonomy to each node using QIIME
# Aligned according to SILVA 16s (since our data is 16s)
align_seqs.py -i data/Westbeach/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras.fasta -t example/path/../../../taxonomy_databases/SILVA_128/16s/rep_set_aligned/99_SILVA_128_aligned_rep_set.fasta -o aligned_seqs

# Filter alignments. Since this is not clade specific, we will include a 5% entropic threshold
filter_alignment.py -i aligned_seqs/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_aligned.fasta -o filter_alignment_G90_E05 -g 0.90 -e 0.05

# Making a phylogeny from alignments (using fasttree)
make_phylogeny.py -i filter_alignment_G90_E05/NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_aligned_pfiltered.fasta -o my_phylogeny.fasttree.tre -t fasttree
# Can check output to see if the tree makes sense

######################
# DIVERSITY ANALYSES #
######################

# ALPHA #
# Rarefaction (note: mapping file is the metadata file)
alpha_rarefaction.py --otu_table_fp otu_table/OTU_Table.wtaxa.biom --output_dir alphadiversity/ --tree_fp my_phylogeny.fasttree.tre --mapping_fp data/Westbeach_Metadata.txt
# The rarefaction plot can be used to assess the appropriate level to rarefy at
# Data must always be rarefied for alpha and beta diversity analyses

# Calculate alpha diversity
alpha_diversity.py -i OTU_table/OTU_table.wtaxa.biom -m chao1,PD_whole_tree -o alphadiversity/adiv_chao1_pd.txt -t my_phylogeny.fasttree.tre

# BETA #
# Create new text file, and write (this is our parameters file):
beta_diversity:metrics unweighted_unifrac,weighted_unifrac,bray_curtis
# Save under data/Westbeach/params.beta.txt

# Calculate Beta diversity
beta_diversity_through_plots.py --otu_table_fp otu_table/OTU_Table.wtaxa.biom --output_dir betadiversity/ --tree_fp my_phylogeny.fasttree.tre --mapping_fp data/Westbeach_Metadata.txt -p data/Westbeach/params.beta.txt 

# Summarize taxa by plots 
summarize_taxa_through_plots.py -i OTU_table/OTU_table.wtaxa.biom -o ./summarize_taxa



























