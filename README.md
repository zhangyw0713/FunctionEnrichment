# FunctionEnrichment
The scripts for GO, KEGG, and MsigDB hallmark enrichment analyses in batch mode.

The specific workflow is shown in the below.

1. For GO and KEGG analyses.

    a. Got the relationships between genes and GO/KEGG terms via the Rscript go.update.R and KEGG.update.R.
    
    b. Download the Gene_info data from NCBI (ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/)
    
    c. Run go.enrich.pl or kegg.enrich.pl to obtain the enriched GO/KEGG IDs.
    
    d. Run get.term.pl to map the IDs to terms and descriptions. 
    
2. For MsigDB analysis.
    
    a. Download the hallmarks from MsigDB database
   
    b. Run MsigDB.R
    
 
