# metagenomic_orobanches
Metagenomic project - in collaboration with Dr. Lucie POULIN and PhD Lisa MARTINEZ (LBPV lab, Nantes, FRANCE)

Here are the codes that I used to make some alpha and beta diversity analyses.

# Summary of the biological context :
Orobanches are holoparasitic plants that are not able to make photosynthesis. These plants are parasitic to many crops (e.g., rapeseeds). 
A farmer has a field that is divided by a road and interestingly, orobanches are mostly growing in one side compared to the other one. 
We want to understand what are the factors of this observation. 
In the lab, they performed some experiments (co-cultures).  

In the code, you can find these abbreviations : 
CN = conducive (sensitive soil) *versus* SN = suppressive (resistant) 


*Need to read it in this way (R script) :*
# First part : 
Do some t-tests on metadata to understand if soils (2 in our case) had an impact on the bacterial community and led to the growth or not of orobanches. 

# Second part :
Alpha diversity (Shannon indices using the phyloseq package) from the abundance table and using the taxa table (not provided here).

# Third part : 
Transform (Centered-Log Ratio) and rarefy the abundance table (also called OTU table) --> two ways to study abundance tables.  

# Fourth part : 
Beta diveristy : 
A. PCoA from a phyloseq object and rarefaction (end of the R script) 
B. UMAP analysis after CLR transformation (python script)
                
Since it is not published yet, I'm only allowed to share codes, not data.  
             
