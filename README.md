# Tritrophic-Interactions
Scripts behind the figures for our paper on the role of a knockout in a key receptor that mediates defenses against herbivory in Phaseolus vulgaris.

In the interest of transparency, these scripts were created in order to replicate the fantastic visual work of Figure 2 from Espina et al., 2024 (10.1111/tpj.17026)
These are by no means perfect, but it accomplished a similar goal.

The way I accomplished this was by using a Bash script to use a BWA to process the fastqc.zip files I received from the service that performed our WGS (Genewiz/Azenta). 
This was aligned to reference genome Labor Ovalle, as our beans are Mesoamerican ("Phaseolus vulgaris Labor Ovalle v1.1 DOE-JGI, https://phytozome-next.jgi.doe.gov/info/PvulgarisLaborOvalle_v1_1").
These scripts were written because we crossed Parent B and a selection of subsequent progeny into Parent A six times, and we wanted to see how big the introgression was. The gene of interest that we introgressed had a deletion in a key receptor.
The idea was that each of our sibling lines should be mostly identical to each other exept for the introgression. 

After the script compiled and deduplicated the BAMs, a subscript kept only primary aligned reads, excluding multi-mapping. 
I ran another script that called single-nucleotide polymorphisms (SNPs) on all chromosomes and compiled them into VCF files. Chromosome 7 carried our gene of interest. 

My third script in R allowed me to turn those SNP calls into a dataframe, where I assigned the SNP calls A-like, B-like, A/B, or Neither. This restricted the plot to informative SNPs only, where the parents were homozygous for opposite SNPs.
Additionally, I wanted to plot SNP density along chromosome 7 as well as depth of coverage for our specific gene of interest, since the SNP calls would not show the deletion.

I was then able to plot those values, export them to PDF, and arrange them into figure panels in Adobe Illustrator. 

Finally, I used a Bash script to compile statistics from all of the BAMs (the initial unfiltered BAMs). An R script then parsed those statistics and compiled them into a summary Excel Spreadsheet. 

The order of the Scripts is 1) the BWA from FASTQ Bash script, 2) The SNP-caller and VCF creator Bash Script, 3) The VCF Analysis in R, 4) the WGS statistic compiler in Bash, 5) the Statistics parser in R.

This was our supplemental figure for our paper documenting the story of how some legumes use this key receptor to detect caterpillar herbivory and mediate direct (chemical) and indirect (volatile for predator) defenses against these herbivores. 
These introgressed lines were our main tool to study this interaction, and this figure gives us a look under the hood of our tools. You can find our preprint here: https://doi.org/10.1101/2025.07.29.667524
