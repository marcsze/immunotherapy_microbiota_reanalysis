REFS = data/references
FIGS = results/figures
TABLES = data/process/tables
PROC = data/process
FINAL = submission/

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'


################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy. Note that this code
# assumes that mothur is in your PATH. If not (e.g. it's in code/mothur/, you
# will need to replace `mothur` with `code/mothur/mothur` throughout the
# following code.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds
	mv trainset14_032015.pds/* $(REFS)/
	rm -rf trainset14_032015.pds
	rm Trainset14_032015.pds.tgz

################################################################################
#
# Part 2: Get and run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################

# Run the mothur analysis to process fastq into a relevent files
# shared, fasta, dist, summary, taxonomy

$(PROC)/kitbias.% : code/Mia.batch code/stabilityedit.R
	bash code/Mia.batch

# Run alpha comparisons

$(TABLES)/kruskal_alpha_pres.csv\
$(TABLES)/kruskal_alpha_ext_kit.csv : $(PROC)/kitbias.groups.ave-std.summary\
code/bias_project_stats.R
	R -e "source('code/bias_project_stats.R')"


# Run alpha relative diffs comparisons

$(TABLES)/kruskal_diffs_w_ExtKit_%\
$(TABLES)/diff_table_w_ExtKit_%\
$(TABLES)/wilcox_diffs_w_Preserve_%\
$(TABLES)/diff_table_w_Preserve_% : $(PROC)/kitbias.groups.ave-std.summary\
$(PROC)/sample_metadata.csv code/run_alpha_rel_diffs_comparison_pres.R\
code/run_alpha_rel_diffs_comparison_ext.R
	R -e "source('code/run_alpha_rel_diffs_comparison_pres.R')"
	R -e "source('code/run_alpha_rel_diffs_comparison_ext.R')"


# Run bray-curits based beta diversity analysis

$(TABLES)/bc_ext_kruskal.csv $(TABLES)/bc_pres_kruskal.csv\
$(TABLES)/bc_ext_dunns.csv $(TABLES)/bc_pres_dunns.csv\
$(TABLES)/bc_raw_graph_data.csv $(TABLES)/bc_summary_graph_data.csv\
$(TABLES)/bc_all_raw_tidy_graph_data.csv $(TABLES)/bc_permanova_analysis.csv\
$(TABLES)/bc_betadisper_analysis.csv $(TABLES)/bc_raw_nmds_values_for_graph.csv\
$(TABLES)/bc_permanova_sample_id_analysis.csv : $(PROC)/kitbias.braycurtis.0.03.lt.ave.dist\
code/run_bray_curtis_preservation_test.R $(PROC)/sample_metadata.csv\
code/run_nmds_permanova_analysis.R
	R -e "source('code/run_bray_curtis_preservation_test.R')"
	R -e "source('code/run_nmds_permanova_analysis.R')"


# Set up specific Taxa variables for preservation 

TAXA = phyla_data family_data genus_data otu_data
TAXA_KRUSK_PRES=$(foreach S, $(TAXA_USED), $(PROC)/$(S)_kruskal.csv)
TAXA_DUNN_PRES=$(foreach S, $(TAXA_USED), $(PROC)/$(S)_dunns.csv)
TAXA_KRUSK_DIFF_PRES=$(foreach S, $(TAXA_USED), $(PROC)/kruskal_diffs_w_ExtKit_$(S).csv)
TAXA_DIFF_PRES=$(foreach S, $(TAXA_USED), $(PROC)/diff_table_w_ExtKit_$(S).csv)


# Run taxa based preservation method analysis

$(TAXA_KRUSK_PRES) $(TAXA_DUNN_PRES) $(TAXA_KRUSK_DIFF_PRES)\
$(TAXA_DIFF_PRES) : $(PROC)/kitbias.0.03.subsample.shared $(PROC)/kitbias.taxonomy\
code/run_taxa_preservation_analysis.R code/run_taxa_rel_diffs_comparison.R
	R -e "source('code/run_taxa_preservation_analysis.R')"
	R -e "source('code/run_taxa_rel_diffs_comparison.R')"


# Set up specific Taxa variables for extraction kit

TAXA_KRUSK_EXT=$(foreach S, $(TAXA_USED), $(PROC)/$(S)_ext_kit_kruskal.csv)
TAXA_DUNN_EXT=$(foreach S, $(TAXA_USED), $(PROC)/$(S)_ext_kit_dunns.csv)
TAXA_GRAPH_EXT=$(foreach S, $(TAXA_USED), $(PROC)/rel_abund_$(S)_graph_data.csv)
TAXA_WILCOX_EXT_PRES=$(foreach S, $(TAXA_USED), $(PROC)/wilcox_diffs_w_Preserve_$(S).csv)
TAXA_DIFF_EXT=$(foreach S, $(TAXA_USED), $(PROC)/diff_table_w_Preserve_$(S).csv)


# Run taxa based extraction kit analysis

$(TAXA_KRUSK_EXT) $(TAXA_DUNN_EXT) $(TAXA_WILCOX_EXT_PRES)\
$(TAXA_DIFF_EXT) $(TAXA_GRAPH_EXT) : $(PROC)/kitbias.0.03.subsample.shared\
$(PROC)/kitbias.taxonomy code/run_taxa_ext_kit_analysis.R\
code/run_taxa_rel_diffs_comparison_ext.R
	R -e "source('code/run_taxa_ext_kit_analysis.R')"
	R -e "source('code/run_taxa_rel_diffs_comparison_ext.R')"


# Set up specific Taxa variables for formed/unformed stool

TAXA_KRUSK_FORM=$(foreach S, $(TAXA_USED), $(PROC)/formed_stool_$(S)_wilcox.csv)


# Run taxa based extraction kit analysis

$(TAXA_KRUSK_FORM) : $(PROC)/kitbias.0.03.subsample.shared $(PROC)/kitbias.taxonomy\
code/run_taxa_formed_stool_analysis.R
	R -e "source('code/run_taxa_formed_stool_analysis.R')"



################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################


#Figure 1. 
# Differences between various preservation methods across extraction kits 
# tested relative to no preservation for Shannon diversity, evenness, and richness.

#Figure 2. 
# Differences between extraction kits across preservation method tested relative to 
# Power Soil extraction kit for Shannon diversity, evenness, and richness.

$(FIGS)/Figure1.pdf\
$(FIGS).Figure2.pdf : $(TABLES)/diff_table_w_ExtKit_diversity.csv\
$(TABLES)/diff_table_w_Preserve_diversity.csv code/make_alpha_diff_graph.R
	R -e "source('code/make_alpha_diff_graph.R')"



#Figure 3.
# Heatmap of significantly different taxa between preservation methods. 

$(FIGS)/Figure3.pdf : $(TABLES)/phyla_data_dunns.csv $(TABLES)/family_data_dunns.csv\
$(TABLES)/genus_data_dunns.csv $(TABLES)/otu_data_dunns.csv\
code/make_dunn_taxa_graph.R
	R -e "source('code/make_dunn_taxa_graph.R')"

#Figure 4. 
# Heatmap of significantly different taxa between extraction kits used. 

$(FIGS)/Figure4.pdf : $(TABLES)/phyla_data_ext_kit_dunns.csv\
$(TABLES)/family_data_ext_kit_dunns.csv $(TABLES)/genus_data_ext_kit_dunns.csv\
$(TABLES)/otu_data_ext_kit_dunns.csv code/make_dunn_ext_kit_taxa_heatmap_graph.R
	R -e "source('code/make_dunn_ext_kit_taxa_heatmap_graph.R')"

#Figure 5. 
# An NMDS plot of all samples using the Bray-Curtis distance metric.

#Figure S1. 
# NMDS plots of all possible combinations of preservation method and 
# extraction kit using the Bray-Curtis distance metric.

$(FIGS)/Figure5.pdf\
$(FIGS)/FigureS1.pdf : $(TABLES)/bc_raw_nmds_values_for_graph.csv\
code/make_bc_formed_stool_nmds_graph.R
	R -e "source('code/make_bc_formed_stool_nmds_graph.R')"

#Figure 6. 
# Differences in bacterial community based on preservation method or 
# extraction kit using the Bray-Curtis distance metric.    

$(FIGS)/Figure6.pdf : $(TABLES)/bc_raw_nmds_values_for_graph.csv code/make_bc_graph.R
	R -e "source('code/make_bc_graph.R')"
	





################################################################################
#
# Part 4: Pull it all together
#
# Render the manuscript
#
################################################################################


write.paper : $(FINAL)/manuscript.Rmd $(FINAL)/supplement.Rmd\
$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
$(FIGS)/Figure3.pdf $(FIGS)/Figure4.pdf\
$(FIGS)/Figure5.pdf $(FIGS)/Figure6.pdf\
$(FIGS)/FigureS1.pdf code/Run_render_paper.R
	R -e "source('code/Run_render_paper.R')"