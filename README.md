# Reciprocal-Best-Hit
A program that takes two blastp output files and filters for reciprocal best hits. Removes duplicate "best hits" and writes the final filtered matches to a summary table containing gene names. (via biomart files)

For this program to work efficiently, sort blastp output files by query ID AND e-value. The commands ran are found in bash-commands file in this repo. 


# To Do
Generalize RBH.py : Currently, variables are named by zebrafish and human- since they were the initial input. Generalize these variables. Perhaps based on argparse naming system? 
