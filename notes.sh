############################
#MED & QIIME basic tutorial#
############################

######
#note: this tutorial begins after data has been demultiplexed, trimmed, filtered, and prepared for use with MED. for an overview of this process, please see steps 1-3 of this lab protocol: https://docs.google.com/document/d/1L83TgSxGLkJ3yyNPexzBB1qtrDlyrj1JwWiDEDjXusc/edit?usp=sharing
######

#our sample data consists of 15 16S barcoded mouse fecal samples from three animal research facilities (5 samples from each facility)
#we'd like to see if there are differences in the microbial communities in each facility

###########
#USING MED#
###########
#the first step after preparing our data is to run MED. 
#MED is invoked like this:
decompose data/seqs.trimmed_filtered_250bp.MED.fna -o MED/decompose-M_50 #this takes about 4 minutes on a 3 core machine, but we'll cheat to make it faster by increasing the minimum substantive abundance per node
decompose data/seqs.trimmed_filtered_250bp.MED.fna -M 50 -o MED/decompose-M_50

#MED sorts sequences into groups (nodes) then calculates the entropy of the nodes based on sequence similarity. 
#each node is divided in two to reduce entropy, and the process is repeated until a minimum threshold of entropy is reached. 
#parameters that can be adjusted include:
# *minimum node size (how many sequences in entire node)
# *minimum node core size (how many sequences belong to node center)
# *maximum variation allowed per node (effects entropy threshold)

#when choosing a non-default M, look closely at your data, at the number of reads per each sample after you do your sequence trimming and filtering.
#generally, it's a good idea to make sure that your -M value won't eliminate taxa that are present in samples that have lower than average numbers of reads.

#example: I want to capture taxa that are only represent 1% of my total community population, and may only be present in a handful of samples. My data are variable in their read counts, some have over 50,000 reads but many only have 10,000 or less.
#I look at my summary of reads per sample, and choose an M value that corresponds to ~1% of the reads for one of my samples that has few reads but still has enough data for me to keep.

##############################
#prepare MED output for QIIME#
##############################
#transpose sample/OTU matrix
ls MED/decompose-M_50/MATRIX-COUNT.txt | perl scripts/transpose.pl -
#MED adds a unique identifier for each node, but leaves leading zeroes in the IDs, so they are all the same character length (ex. 001, 058), we remove leading zeros before using QIIME
sed 's/^0*//' MED/decompose-M_50/MATRIX-COUNT.txt_transposed > MED/decompose-M_50/MATRIX-COUNT.transposed.txt

#MED also includes leading zeros in it's node representative fasta file, and the number of reads that belong to each node. QIIME doesn't like either of these things, so we remove them
sed 's/|size:[0-9]*$//g' MED/decompose-M_50/NODE-REPRESENTATIVES.fasta | sed 's/^>0*/>/' > MED/decompose-M_50/NODE-REPRESENTATIVES.DOWNSTREAM.fasta

#To install qiime for macs:
macqiime

####################
#ASSIGNING TAXONOMY#
####################
#using NODE-REPRESENTATIVES.DOWNSTREAM.fasta to assign taxonomy to each node using QIIME
assign_taxonomy.py -i MED/decompose-M_50/NODE-REPRESENTATIVES.DOWNSTREAM.fasta -t db_files/88_SILVA_128_taxa_map_7_levels.txt -r db_files/88_SILVA_128_prks_rep_set.fasta -m uclust -o ./assign_taxonomy

##############################
#MAKING THE PHYLOGENETIC TREE#
##############################
#now align sequences #takes about 2 min
align_seqs.py -i MED/decompose-M_50/NODE-REPRESENTATIVES.DOWNSTREAM.fasta -t db_files/88_SILVA_128_aligned_rep_set.fasta -o aligned_seqs

#now filter the alignment #remove positions with >90% gaps and the 5% most entropic positions
#-g is a gap filter and -e is an entropy filter. -e also suppresses the lane mask (which is enabled by default). the lane mask is for the greengenes alignment (as far as I am aware), so if you don't want to filter by entropy, you must supress it with -s.
filter_alignment.py -i aligned_seqs/NODE-REPRESENTATIVES.DOWNSTREAM_aligned.fasta -o filter_alignment_G90_E05 -g 0.90 -e 0.05

#now make the phylogenetic tree
# -t is the method for tree building. I use fasttree for the basic analysis
make_phylogeny.py -i filter_alignment_G90_E05/NODE-REPRESENTATIVES.DOWNSTREAM_aligned_pfiltered.fasta -o 16s_makephylo_fasttree.tre -t fasttree

######################
#MAKING THE OTU TABLE#
######################
mkdir OTU_table
biom convert -i MED/decompose-M_50/MATRIX-COUNT.transposed.txt -o OTU_table/OTU_Table.biom --to-json --table-type="OTU table"

#summarize it
biom summarize-table -i OTU_table/OTU_Table.biom > OTU_table/OTU_Table.summary
#viewing the summary file allows us to check and see if we need to remove any samples for lack of reads, and to make sure our negative controls were successful

#we can create parameters files that control what calculations QIIME does during alpha and beta diversity analyses
echo "alpha_diversity.py:metrics chao1,PD_whole_tree,equitability" > alpha_params.txt
echo "beta_diversity:metrics bray_curtis,unweighted_unifrac,weighted_unifrac" > beta_params.txt

#alpha diversity
alpha_rarefaction.py --otu_table_fp OTU_table/OTU_Table.biom --output_dir alphadiversity/ --tree_fp 16s_makephylo_fasttree.tre --mapping_fp metadata/mouse_cortisol.mapping_file.txt -p alpha_params.txt
#based on alpha rarefaction plots we can remove samples with fewer than 8301 reads (none of them)

#now we add the taxonomy
biom add-metadata -i OTU_table/OTU_Table.biom -o OTU_table/OTU_Table.wtaxa.biom --observation-metadata-fp assign_taxonomy/*.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy
cd .. #back to main directory

#calculate beta diversity
beta_diversity_through_plots.py --otu_table_fp OTU_table/OTU_Table.wtaxa.biom --output_dir betadiversity/ --tree_fp 16s_makephylo_fasttree.tre --mapping_fp metadata/mouse_cortisol.mapping_file.txt -p beta_params.txt

#summarize taxa in OTU table
summarize_taxa_through_plots.py -i OTU_table/OTU_Table.wtaxa.biom -o ./summarize_taxa
