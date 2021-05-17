## 2021/05/17
In this point we are going to archive the project until further notice.

### RNAscope
We are going to apply RNAscope of Smim41 and Col9a2 to detect the kranocyte populations, or whichever population that is expressing those markers.
At first I tried to select the best region for RNAscope by aligning the reads and selecting the ones mapping to the genes, but then talking to Olga I realized that the reads were 3' biased and did not have UMI deduplication, so it 
was highly unlikely that the region of interest was going to be there. Also, we found that the RNAscope page had an app to build the probes, so it will be easier fot th ecompany to do it.

### Rest of the project
In the Giordani and Dell'Orso datasets we found that, in FAPs, there were two krano A populations: one in an smal *independent* cluster, and a second one *inside* the FAP cluster (info in 1_ notebook). These populations are not reproducible in the rest of the datasets, but being reproducible in two, it is an interesting finding. Once the RNAscope is done I think we should go further into these populations. 