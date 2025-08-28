# Tritrophic-Interactions
Scripts behind the figures for our paper on the role of a knockout in a key receptor that mediates defenses against herbivory in Phaseolus vulgaris.

In the interest of transparency, these scripts were created in order to replicate the fantastic visual work of Figure 2 from Espina et al., 2024 (10.1111/tpj.17026)
These are by no means perfect, but it accomplished the same goal. 

The way I accomplished this was by using a Bash script to use a BWA to process the fastqc.zip files I received from the service that performed our WGS (Genewiz/Azenta). This was aligned to reference genome G19833. 
Even though our beans are Mesoamerican, we chose to align it to this Andean version because all that mattered was whether the SNPs were more like Parent A or Parent B, and more work has been done using G19833 than with the Mesoamerican refernce Labor Ovalle. 
This was because we introgressed Parent B into Parent A six times, and we wanted to see how big the introgression was. The gene of interest that we introgressed had a deletion in a key receptor.
The idea was that each of our sibling lines should be identical to each other exept for the introgression. 

After the script compiled and deduplicated the BAMs, I ran another script that specifically called SNPs from Chromosome 7, which carried our gene of interest. 
This script compiled the calls of all 8 lines, the 2 parents and the 3 sibling lines. This gave me a key output - a VCF file. 

My third and final script in R allowed me to turn those SNP calls into a dataframe, where I assigned the SNP calls A-like, B-like, A/B, or Neither. 
Additionally, I wanted to plot SNP density along chromosome 7 as well as depth of coverage for our specific gene of interest, since the SNP calls would not show the deletion. 

I was then able to plot those values, export them to EPS, and arrange them into figure panels in Illustrator. 
