 The scripts analyse Nanopore sequences from HIV drug-resistance specimens. The scripts are run in the order they are numbered, and in some cases the output of one script becomes the input to the next. All input files are included. Note should be taken to edit the input directory where necessary for each script to match the user path to the input file directory. 

1.	Dorado basecalls bash script on the HPC (Slurm using GPU). 
Usage:	create the bash script using the nano command (script in folder scripts)
	nano 01_dorado_basecall.slurm 
	Exit the bash file: ctrl-x
	Save the bash file: y
	Esc		
	Run: sbatch 01_dorado_basecall.slurm	
	Check the progress of your job: squeue -u your_username 

		Script: 01_dorado_basecall_gpu_r10_model.slurm
 
2.	Demultiplex the Dorado basecalled reads using a Slurm script created and run as in 1 above. 
02_dorado_demux_gpu_batch.slurm

3.	Remove human and SARS-CoV-2 genome reads using removehost you can install from the Github repository (https://github.com/jiangweiyao/SanitizeMe)
Note: The script “removehost.py” in the execution bin was modified so that the output files could preserve the timestamp of the fasta file by creating an intermediate sequence ID file and saving it as script 03_removehost.py. It is recommended that you save it with the original name- removehost.py before overwriting the original file.  

		Script: 03_removehost.py 

4.	Filter the reads by length to retain 500-1200 base pairs and quality  score 15 using Nanoq that can be installed from GitHub at https://github.com/esteinig/nanoq. To iterate across all files in the folder, a bash script 04_nanoq.sh is included. Run the script in the terminal by making it executable.

Chmod +x 04_nanoq.sh
./04_nanoq.sh or bash 04_nanoq.sh

5.	Perform additional polishing by running Dechat that can be installed from the Github repository https://github.com/LuoGroup2023/DeChat.

Note: kmers, k, were set to 35, the default is 21 and maximum 64. The Dechat bash script iterates over the folder and runs all, “05_dechat_restore_headers_iterate.slurm”. This was run on a remote high-performance cluster (HPC).

		sbatch 05_dechat_restore_headers_iterate. 

6.	Orient the Sanger generated consensus fasta files as the forward (sense) strand using a Python script “06_ seq_reverse_complement_iterate.py”. The input fasta files are in the folder / fasta_pr and the output forward oriented fastas /fasta_pr_references_fowardclear.

Python slurm06_ seq_reverse_complement_iterate.py


7.	Reference-based assembly of the fastq files to the corresponding Sanger sequences  as references using Minimap2 https://github.com/lh3/minimap2  and a bash script, “07_minimap2_iterate_refences.sh” to iterate through the folders matching a reference to its corresponding Nanopore reads in the folder /output_dechat/barcode*. Make the bash script executable.

chmod +x 07_minimap2_iterate_refences.sh
./07_minimap2_iterate_refences.sh

8.	Polish the output BAM files using Medaka https://github.com/nanoporetech/medaka and a shell script to iterate through the BAM files generated in 7 above and their corresponding Sanger reference sequences. Make the bash script executable.

chmod +x 08_medaka_cpu.sh
./08_medaka_cpu.sh

9.	Repeat the Minimap2 reference-based assembly using the Medaka-generated consensus as the reference using the same script in 7 above, adjusting the path to the references (07_minimap2_iterate_refences.sh).
	

10.	Variant call using Claire3 https://github.com/HKU-BAL/Clair3 . Use the super accuracy model in the folder named /clair3_models. Use a bash script to iterate across BAM files; “09_clair3_iterate_medaka.sh”. Make the script executable.

chmod +x 09_clair3_iterate_medaka.sh
09_clair3_iterate_medaka.sh

11.	Split the BAM files into various sub-samples for re-running Clair3 to determine the level where variants level off in number (addition of sequences does not improve the number of variants). Use the bash script  “10_clair3_iterate4_medaka_split_bams.sh”.
  
chmod +x 10_clair3_iterate4_medaka_split_bams.sh
./10_clair3_iterate4_medaka_split_bams.sh

12.	HIV Drug Resistance: the HyDra was used for the pol region - https://hydra.canada.ca/analyses/new?lang=en-CA. The input parameters were;
consensus per cent 20%
target coverage 5000
length cutoff 1000
score cutoff 15
minimum read depth 1000
minimum variant quality 15
error rate 0.05
minimum allele count 5
Minimum AA Frequency 0.01

13.	Determine HIV drug resistance profiles for both Sanger sequences and matched Nanopore sequences using the Stanford HIV drug resistance database online tool  https://hivdb.stanford.edu/hivdb/by-reads/.

Note: The NGS Nanopore needs the fastq files to be converted into codon frequencies (CodFreq files) that are then input into the database. The input fastq files were generated in step 4 above (those filtered with Nanoq) .

Statistics and Plotting Figures

14.	The sequence depth scatter plot data is in the folder “/depths_minimap2”. Each sample depth per position was derived and stored as a txt file. The files were combined into a single text file, “all_depths.txt” which is used by the Python script “11_depth.py” and creates a scatter plot “average_depth_all_samples_absolute_scale_skipped_labels.png”

python 11_depth.py 

15.	Viral Load Plots: 
Barplot
input file “viral_load.csv”
Script “12_viral_load_barchart.py”

Scatter plots: viral load versus sequence depth and mutation burden
	Input files “sequenceSummaries_nanopore.csv” and “viral_load.csv”	        	 
	Scripts “13_viral_load_depth.” and  “py14_viral_load_mutation.py”.

16.	Comparing Nanopore versus Sanger HIVDR Levels concordance (Bland-Altman plot, McNemar’s test for paired variables, Overall per cent agreement, Heatmap visualisation of concordance, etc). The python script “15_nanopore_vs_sanger_HIVDR2.py” was used. The script requires 2 csv file inputs from the Stanford HIV database output i.e.
1.	resistanceSummaries_nanopore.csv
2.	resistanceSummaries_sanger.csv  

python 15_nanopore_vs_sanger_HIVDR2.py 

17.	Sanger vs Nanopore (Hydraweb) HIVDR concordance performs a position-level comparison of HIV mutation calls derived from Sanger sequencing and Oxford Nanopore sequencing for Stanford database-derived versus Hydraweb HIVDR profiles (Jaccard Similarity). The python script “python 16_compareDR_hydraweb_sanger.py” expects 2 folders:
1.	/mutation-list_nanopore_stanford
2.	/mutation-list_Sanger_stanford

python 16_compareDR_hydraweb_sanger.py 

18.	Sanger vs Nanopore (Stanford) codon concordance was determined using Sanger sequences and the minimap2 assembly results from Nanopore sequences using 2 Python scripts. The 2 files were input into the script “17_codon_bias_from_bam.py “ and generated an output file “global_comparison_summary.csv” that was then input for processing the heatmap into the script “18_codon_freq_Heatmap.py”.
python 17_codon_bias_from_bam.py
python 18_codon_freq_Heatmap.py

19.	Depth distribution of variants across the down-sampled data was analysed using a Bash script “19_split_bam.sh” that splits the BAM files and saves them in a folder “/subsampled_bams”. A Python script “20_compare_subsamples.py” then analyses the subsampled BAM files to visualise how variant-calling performance (F1 Score) changes across different BAM sampling proportions (i.e. 0.005, 0.01, 0.02, 0.05 and 0.1). This generates a log-scaled boxplot. The input BAM files are in the folder “/minimap2_results2_medaka”. 

./19_split_bam.sh
python 20_compare_subsamples.py

20.	Subtype concordance was determined by comparing the Stanford HIVDR database-assigned subtype for the Sanger-derived sequences with the Nanopore-derived sequences. The 2 input files from the Stanford HIVDR output are: “sequenceSummaries_sanger.csv” and “sequenceSummaries_nanopore.csv”. The Python script “21_subtypes_matrix.py” generates a confusion matrix.

python 21_subtypes_matrix.py
