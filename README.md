# HMMerProteinSeq



> User can use protein accession (will be downloaded automatically from NCBI RefSeq database) or protein sequence in a fasta format as input. 
> This script uses HMMs identified by Chayan Kumar Saha and Gemma Atkinson; used in the analyses of Jimmy et al PNAS 2020.
> Also, user can put any HMM model of their own in the HMM directory to scan the query protein sequence. 



Command: 

If input is a protein accession (for example, WP_090558406.1):
>    python3 HMMerProteinSeq.py -a WP_090558406.1 -u userid@email.com -d HMMs/

If input is a protein sequence in fasta format in a text file (for example, WP_090558406.1.fasta or WP_090558406.1.faa):
>    python3 HMMerProteinSeq.py -f WP_090558406.1.fasta -u userid@email.com -d HMMs/
>    or, python3 HMMerProteinSeq.py -f WP_090558406.1.faa -u userid@email.com -d HMMs/

