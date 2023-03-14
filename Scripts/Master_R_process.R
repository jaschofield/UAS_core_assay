## Master R script to process clustered barcode sequences and create final data table
## Jeremy Schofield, Hahn Lab. 2023.

library(tidyverse)
library(purrr)
library(Biostrings)
library(fuzzyjoin)
Date <- Sys.Date()

## Read formatted barcode files output from second clustering script
MED15_plus_1 <- read.table('MED15_plus_clust_table_1.txt')
MED15_plus_2 <- read.table('MED15_plus_clust_table_2.txt')
MED15_plus_3 <- read.table('MED15_plus_clust_table_3.txt')
MED15_plus_4 <- read.table('MED15_plus_clust_table_4.txt')
MED15_minus_1 <- read.table('MED15_minus_clust_table_1.txt')
MED15_minus_2 <- read.table('MED15_minus_clust_table_2.txt')
MED15_minus_3 <- read.table('MED15_minus_clust_table_3.txt')
MED15_minus_4 <- read.table('MED15_minus_clust_table_4.txt')
SPT7_plus_1 <- read.table('SPT7_plus_clust_table_1.txt')
SPT7_plus_2 <- read.table('SPT7_plus_clust_table_2.txt')
SPT7_plus_3 <- read.table('SPT7_plus_clust_table_3.txt')
SPT7_plus_4 <- read.table('SPT7_plus_clust_table_4.txt')
SPT7_minus_1 <- read.table('SPT7_minus_clust_table_1.txt')
SPT7_minus_2 <- read.table('SPT7_minus_clust_table_2.txt')
SPT7_minus_3 <- read.table('SPT7_minus_clust_table_3.txt')
SPT7_minus_4 <- read.table('SPT7_minus_clust_table_4.txt')
TAF13_plus_1 <- read.table('TAF13_plus_clust_table_1.txt')
TAF13_plus_2 <- read.table('TAF13_plus_clust_table_2.txt')
TAF13_plus_3 <- read.table('TAF13_plus_clust_table_3.txt')
TAF13_plus_4 <- read.table('TAF13_plus_clust_table_4.txt')
TAF13_minus_1 <- read.table('TAF13_minus_clust_table_1.txt')
TAF13_minus_2 <- read.table('TAF13_minus_clust_table_2.txt')
TAF13_minus_3 <- read.table('TAF13_minus_clust_table_3.txt')
TAF13_minus_4 <- read.table('TAF13_minus_clust_table_4.txt')

## Create table for each file pulling out barcode, UMI, and cluster sizes

RNA_MED15_1_plus <- MED15_plus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 1)

RNA_MED15_2_plus <- MED15_plus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 2)

RNA_MED15_3_plus <- MED15_plus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 3)

 RNA_MED15_4_plus <- MED15_plus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 4)

RNA_MED15_1_minus <- MED15_minus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 1)

RNA_MED15_2_minus <- MED15_minus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 2)

RNA_MED15_3_minus <- MED15_minus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 3)

 RNA_MED15_4_minus <- MED15_minus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "MED15") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 4)

RNA_SPT7_1_minus <- SPT7_minus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 1)

RNA_SPT7_2_minus <- SPT7_minus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 2)

RNA_SPT7_3_minus <- SPT7_minus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 3)

 RNA_SPT7_4_minus <- SPT7_minus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 4)

RNA_SPT7_1_plus <- SPT7_plus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 1)

RNA_SPT7_2_plus <- SPT7_plus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 2)

RNA_SPT7_3_plus <- SPT7_plus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 3)

 RNA_SPT7_4_plus <- SPT7_plus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "SPT7") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 4)


RNA_TAF13_1_minus <- TAF13_minus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 1)

RNA_TAF13_2_minus <- TAF13_minus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 2)

RNA_TAF13_3_minus <- TAF13_minus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 3)

 RNA_TAF13_4_minus <- TAF13_minus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "no") %>%
  mutate(replicate = 4)

RNA_TAF13_1_plus <- TAF13_plus_1  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 1)

RNA_TAF13_2_plus <- TAF13_plus_2  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 2)

RNA_TAF13_3_plus <- TAF13_plus_3  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 3)

 RNA_TAF13_4_plus <- TAF13_plus_4  %>%
 separate(V1, c("drop1", "size"), sep = ';size=') %>%
  separate(drop1, c("drop2", "seqs"), sep = ';seqs=') %>%
  separate(V2, c("bc1", "dups"), sep = 'CGAAAGA.............TCAAGAC') %>%
  separate(bc1, c("prim1", "barcode"), sep = 'AGCAGG') %>%
  separate(dups, c("pcrdup", "prim2"), sep = 'CTGTCT') %>%
  select(barcode, pcrdup, size, seqs) %>%
  filter(nchar(barcode) == 18) %>%
  mutate(strain = "TAF13") %>%
  mutate(auxin = "yes") %>%
  mutate(replicate = 4)


barcode_master <- bind_rows(RNA_MED15_1_plus, RNA_MED15_2_plus, RNA_MED15_3_plus, RNA_MED15_4_plus, RNA_MED15_1_minus, 
	RNA_MED15_2_minus, RNA_MED15_3_minus, RNA_MED15_4_minus, RNA_SPT7_1_plus, RNA_SPT7_2_plus, RNA_SPT7_3_plus, RNA_SPT7_4_plus, 
	RNA_SPT7_1_minus, RNA_SPT7_2_minus, RNA_SPT7_3_minus, RNA_SPT7_4_minus, RNA_TAF13_1_plus, RNA_TAF13_2_plus, RNA_TAF13_3_plus, 
	RNA_TAF13_4_plus, RNA_TAF13_1_minus, RNA_TAF13_2_minus, RNA_TAF13_3_minus, RNA_TAF13_4_minus)

saveRDS(barcode_master, paste0("barcode_master_", Date, ".rds"))

## Read formatted plasmid DNA files output from second clustering script

MED_DNA_1 <- read.table('MED15_clust_table_1.txt')
MED_DNA_2 <- read.table('MED15_clust_table_2.txt')
MED_DNA_3 <- read.table('MED15_clust_table_3.txt')
MED_DNA_4 <- read.table('MED15_clust_table_4.txt')
SPT_DNA_1 <- read.table('SPT7_clust_table_1.txt')
SPT_DNA_2 <- read.table('SPT7_clust_table_2.txt')
SPT_DNA_3 <- read.table('SPT7_clust_table_3.txt')
SPT_DNA_4 <- read.table('SPT7_clust_table_4.txt')
TAF_DNA_1 <- read.table('TAF13_clust_table_1.txt')
TAF_DNA_2 <- read.table('TAF13_clust_table_2.txt')
TAF_DNA_3 <- read.table('TAF13_clust_table_3.txt')
TAF_DNA_4 <- read.table('TAF13_clust_table_4.txt')

## Create table for each file pulling out barcode, UAS region, and core promoter region

MED_DNA_1_processed <- MED_DNA_1 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 1, strain = "MED15")

MED_DNA_2_processed <- MED_DNA_2 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 2, strain = "MED15")

MED_DNA_3_processed <- MED_DNA_3 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 3, strain = "MED15")

MED_DNA_4_processed <- MED_DNA_4 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 4, strain = "MED15")

SPT_DNA_1_processed <- SPT_DNA_1 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 1, strain = "SPT7")

SPT_DNA_2_processed <- SPT_DNA_2 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 2, strain = "SPT7")

SPT_DNA_3_processed <- SPT_DNA_3 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 3, strain = "SPT7")

SPT_DNA_4_processed <- SPT_DNA_4 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 4, strain = "SPT7")

TAF_DNA_1_processed <- TAF_DNA_1 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 1, strain = "TAF13")

TAF_DNA_2_processed <- TAF_DNA_2 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 2, strain = "TAF13")

TAF_DNA_3_processed <- TAF_DNA_3 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 3, strain = "TAF13")

TAF_DNA_4_processed <- TAF_DNA_4 %>%
    mutate(UAS_seq = substring(V2, 32, 83), BC_reg = substring(V2, 125, 161), core_reg = substring(V2, 165, 200)) %>%
    separate(BC_reg, c("drop", "BC1"), sep = 'TTCG') %>%
    separate(BC1, c("BC_rev", "drop2"), sep = 'CCTG') %>%
    separate(core_reg, c("drop3", "core_rev"), sep = 'GTCTTG') %>%
    mutate(replicate = 4, strain = "TAF13")


DNA_master <- bind_rows(MED_DNA_1_processed, MED_DNA_2_processed, MED_DNA_3_processed, MED_DNA_4_processed,
    SPT_DNA_1_processed, SPT_DNA_2_processed, SPT_DNA_3_processed, SPT_DNA_4_processed,
    TAF_DNA_1_processed, TAF_DNA_2_processed, TAF_DNA_3_processed, TAF_DNA_4_processed) %>%
  filter(!is.na(BC_rev)) %>%
  filter(!is.na(core_rev)) %>%
    filter(BC_rev > 0 & core_rev > 0)

## Format barcode and core regions as character strings and reduce table size

DNA_master$BC <- reverseComplement(DNAStringSet(DNA_master$BC_rev))
DNA_master$barcode <- as.character(DNA_master$BC)
DNA_master$core <- reverseComplement(DNAStringSet(DNA_master$core_rev))
DNA_master$core_seq <- as.character(DNA_master$core)

DNA_master <- DNA_master %>%
  separate(V1, c("drop_col", "numbers"), sep = ';size=') %>%
  separate(numbers, c("size", "seqs"), sep = ';seqs=') %>%
  dplyr::select(UAS_seq, core_seq, replicate, strain, barcode, size, seqs)

saveRDS(DNA_master, paste0("DNA_master_", Date, ".rds"))

## create joining table with core promoter names and identifying sequences

core_join_ID <- c("KRS1", "NOP13", "RPC10", "SIT1")
core_seq <- c("AATCCGTAATAAACTTCAATAGCA", "CAACTTTCTATAATTTGAGATAGA", "GAATAAAAGTCCACAAGTATAACA", "GACCCAGTACGAAAATTTTCCATA")
core_df <- data.frame(core_join_ID, core_seq)

DNA_master_filt <- DNA_master %>% filter(nchar(barcode) == 18)

## join DNA master file with core promoter table to identify core for each DNA sequence
DNA_master_corejoined <- fuzzyjoin::stringdist_left_join(DNA_master_filt, core_df, by = "core_seq", distance_col = "core_distance", max_dist = 8) %>% mutate(sequence_ID = row_number())

## bring in sequences of UASs from Twist oligo order and create identifying sequence for join
reference_table <- read.csv(file = '//fh/fast/hahn_s/user/jschofie/2206_UAS_redo.dir/UAS_reference_table.csv')
reference_table <- reference_table %>% mutate(UAS_seq = substring(V2, 3, 54))

n <- 20000
nr <- nrow(DNA_master_corejoined)
splits <- split(DNA_master_corejoined, rep(1:ceiling(nr/n), each=n, length.out=nr))

datalist = list()

for (i in 1:ceiling(nr/n)) {
  join <- stringdist_left_join(splits[[i]], reference_table, max_dist = 15, distance_col = "distance")
  datalist[[i]] <- join
  print(i)
}

DNA_UAS_master <- bind_rows(datalist)

saveRDS(DNA_UAS_master, paste0("DNA_master_", Date, ".rds"))


### fill in non-matching core observations, make distance longer than true matches (20)
DNA_UAS_master <- DNA_UAS_master %>%
  mutate_at(8, ~replace_na(.,'none')) %>%
  mutate_at(13, ~replace_na(.,'none')) %>%
  mutate_at(14, ~replace_na(.,'none')) %>%
  mutate_at(15, ~replace_na(.,'none')) %>%
  mutate_at(17, ~replace_na(.,20)) %>%
  mutate(replicate_batch = ifelse(replicate %in% c(1, 2), "A", "B"))

DNA_UAS_master_filt <- DNA_UAS_master %>% group_by(sequence_ID) %>% filter(distance == min(distance)) %>% 
  separate(V1, c("UAS_name", "seqinfo"), sep = "::") %>% mutate(strand = ifelse(grepl('(+)', seqinfo, fixed = TRUE), "plus",
                                                                             ifelse(grepl('(-)', seqinfo, fixed = TRUE), "minus", NA)))


DNA_UAS_master_filt$DNA_size <- as.numeric(DNA_UAS_master_filt$size)
UAS_min_filt_forjoin <- DNA_UAS_master_filt %>%
  group_by(strain, core_join_ID, UAS_name, barcode, replicate_batch) %>% summarize(DNA_ct = mean(DNA_size)) %>% ungroup() %>% 
  dplyr::select(strain, core_join_ID, UAS_name, DNA_ct, barcode, replicate_batch)

colnames(UAS_min_filt_forjoin) <- c("strain", "core_ID", "UAS_name", "DNA_size", "barcode", "replicate_batch")

saveRDS(UAS_min_filt_forjoin, paste0("UAS_filt_master_", Date, ".rds"))

## join barcode master with filtered UAS master
barcode_master$RNA_size <- as.numeric(barcode_master$size)

barcode_master_counts <- barcode_master %>%
  group_by(strain, auxin, replicate, barcode) %>%
  summarize(barcode_count = sum(RNA_size))

UAS_count_table <- full_join(barcode_master_counts, UAS_min_filt_forjoin, by = c("barcode", "strain"))

saveRDS(UAS_count_table, paste0("count_table_", Date, ".rds"))

## perform DNA counts, filter for most common DNA observation, and remove any barcodes appearing in more than one UAS+core context
UAS_counted <- UAS_min_filt_forjoin %>% group_by(strain, core_ID, UAS_name, replicate_batch) %>% summarize(DNA_ct = sum(DNA_size))

UAS_count_filt <- UAS_count_table %>% ungroup() %>% group_by(strain, auxin, replicate, barcode) %>% filter(DNA_size == max(DNA_size))
barcode_multi_list <- UAS_count_filt %>% ungroup() %>% group_by(UAS_name, core_ID, barcode) %>% summarize() %>% 
  group_by(barcode) %>% summarize(n = n())
barcode_remove <- barcode_multi_list %>% filter(n>1)
UAS_count_bc_filt <- UAS_count_filt %>% ungroup() %>% filter(!barcode %in% barcode_remove$barcode)
UAS_count_collapse <- UAS_count_bc_filt %>% ungroup() %>%
  group_by(strain, auxin, replicate, core_ID, UAS_name) %>% summarize(RNA_ct = sum(barcode_count), DNA_ct = sum(DNA_size))

UAS_counted_join <- UAS_count_collapse %>% dplyr::select(-DNA_ct) %>% mutate(replicate_batch = ifelse(replicate %in% c(1, 2), "A", "B"))

## join final barcode counts with DNA counts
UAS_counted_join <- left_join(UAS_counted_join, UAS_counted)

UAS_count_wide <- UAS_counted_join %>% dplyr::select(-replicate_batch) %>% pivot_wider(names_from = c("auxin", "strain", "replicate"), values_from = c("RNA_ct", "DNA_ct")) %>% dplyr::select(-contains("RNA_ct_NA"))

## normalize barcode counts by library size
UAS_count_norm <- UAS_count_wide %>% mutate(RNA_ct_yes_TAF13_1 = RNA_ct_yes_TAF13_1 * 2.672105935,
                                            RNA_ct_yes_TAF13_2 = RNA_ct_yes_TAF13_2 * 2.009443024,
                                            RNA_ct_yes_TAF13_3 = RNA_ct_yes_TAF13_3 * 1.216940713,
                                            RNA_ct_yes_TAF13_4 = RNA_ct_yes_TAF13_4 * 1.254215812,
                                            RNA_ct_yes_SPT7_1 = RNA_ct_yes_SPT7_1 * 6.187131763,
                                            RNA_ct_yes_SPT7_2 = RNA_ct_yes_SPT7_2 * 2.844058115,
                                            RNA_ct_yes_SPT7_3 = RNA_ct_yes_SPT7_3 * 4.064769696,
                                            RNA_ct_yes_SPT7_4 = RNA_ct_yes_SPT7_4 * 2.742724737,
                                            RNA_ct_yes_MED15_1 = RNA_ct_yes_MED15_1 * 2.064694576,
                                            RNA_ct_yes_MED15_2 = RNA_ct_yes_MED15_2 * 2.734977904,
                                            RNA_ct_yes_MED15_3 = RNA_ct_yes_MED15_3 * 2.148153329,
                                            RNA_ct_yes_MED15_4 = RNA_ct_yes_MED15_4 * 1.564525844,
                                            RNA_ct_no_TAF13_1 = RNA_ct_no_TAF13_1 * 2.682537938,
                                            RNA_ct_no_TAF13_2 = RNA_ct_no_TAF13_2 * 1,
                                            RNA_ct_no_TAF13_3 = RNA_ct_no_TAF13_3 * 1.499447138,
                                            RNA_ct_no_TAF13_4 = RNA_ct_no_TAF13_4 * 1.22907855,
                                            RNA_ct_no_SPT7_1 = RNA_ct_no_SPT7_1 * 2.758886711,
                                            RNA_ct_no_SPT7_2 = RNA_ct_no_SPT7_2 * 4.290922989,
                                            RNA_ct_no_SPT7_3 = RNA_ct_no_SPT7_3 * 2.23735417,
                                            RNA_ct_no_SPT7_4 = RNA_ct_no_SPT7_4 * 2.969323077,
                                            RNA_ct_no_MED15_1 = RNA_ct_no_MED15_1 * 2.056229772,
                                            RNA_ct_no_MED15_2 = RNA_ct_no_MED15_2 * 2.177930024,
                                            RNA_ct_no_MED15_3 = RNA_ct_no_MED15_3 * 1.519607569,
                                            RNA_ct_no_MED15_4 = RNA_ct_no_MED15_4 * 1.577304321,
                                            DNA_ct_no_TAF13_1 = DNA_ct_no_TAF13_1 * 1.244320483,
                                            DNA_ct_yes_TAF13_1 = DNA_ct_yes_TAF13_1 * 1.244320483,
                                            DNA_ct_no_TAF13_2 = DNA_ct_no_TAF13_2 * 1,
                                            DNA_ct_yes_TAF13_2 = DNA_ct_yes_TAF13_2 * 1,
                                            DNA_ct_no_TAF13_3 = DNA_ct_no_TAF13_3 * 1.076844887,
                                            DNA_ct_yes_TAF13_3 = DNA_ct_yes_TAF13_3 * 1.076844887,
                                            DNA_ct_no_TAF13_4 = DNA_ct_no_TAF13_4 * 1.32500755,
                                            DNA_ct_yes_TAF13_4 = DNA_ct_yes_TAF13_4 * 1.32500755,
                                            DNA_ct_no_SPT7_1 = DNA_ct_no_SPT7_1 * 1.000396823,
                                            DNA_ct_yes_SPT7_1 = DNA_ct_yes_SPT7_1 * 1.000396823,
                                            DNA_ct_no_SPT7_2 = DNA_ct_no_SPT7_2 * 1.421963579,
                                            DNA_ct_yes_SPT7_2 = DNA_ct_yes_SPT7_2 * 1.421963579,
                                            DNA_ct_no_SPT7_3 = DNA_ct_no_SPT7_3 * 1.200118117,
                                            DNA_ct_yes_SPT7_3 = DNA_ct_yes_SPT7_3 * 1.200118117,
                                            DNA_ct_no_SPT7_4 = DNA_ct_no_SPT7_4 * 1.110036472,
                                            DNA_ct_yes_SPT7_4 = DNA_ct_yes_SPT7_4 * 1.110036472,
                                            DNA_ct_no_MED15_1 = DNA_ct_no_MED15_1 * 1.091716072,
                                            DNA_ct_yes_MED15_1 = DNA_ct_yes_MED15_1 * 1.091716072,
                                            DNA_ct_no_MED15_2 = DNA_ct_no_MED15_2 * 1.103233865,
                                            DNA_ct_yes_MED15_2 = DNA_ct_yes_MED15_2 * 1.103233865,
                                            DNA_ct_no_MED15_3 = DNA_ct_no_MED15_3 * 1.285770184,
                                            DNA_ct_yes_MED15_3 = DNA_ct_yes_MED15_3 * 1.285770184,
                                            DNA_ct_no_MED15_4 = DNA_ct_no_MED15_4 * 1.093085503,
                                            DNA_ct_yes_MED15_4 = DNA_ct_yes_MED15_4 * 1.093085503)

## normalize barcode counts by DNA input
UAS_count_norm <- UAS_count_norm %>% ungroup() %>% mutate(TAF13_plus_1 = RNA_ct_yes_TAF13_1 / DNA_ct_yes_TAF13_1,
                                                              TAF13_plus_2 = RNA_ct_yes_TAF13_2 / DNA_ct_yes_TAF13_2,
                                                              TAF13_plus_3 = RNA_ct_yes_TAF13_3 / DNA_ct_yes_TAF13_3,
                                                              TAF13_plus_4 = RNA_ct_yes_TAF13_4 / DNA_ct_yes_TAF13_4,
                                                              SPT7_plus_1 = RNA_ct_yes_SPT7_1 / DNA_ct_yes_SPT7_1,
                                                              SPT7_plus_2 = RNA_ct_yes_SPT7_2 / DNA_ct_yes_SPT7_2,
                                                              SPT7_plus_3 = RNA_ct_yes_SPT7_3 / DNA_ct_yes_SPT7_3,
                                                              SPT7_plus_4 = RNA_ct_yes_SPT7_4 / DNA_ct_yes_SPT7_4,
                                                              MED15_plus_1 = RNA_ct_yes_MED15_1 / DNA_ct_yes_MED15_1,
                                                              MED15_plus_2 = RNA_ct_yes_MED15_2 / DNA_ct_yes_MED15_2,
                                                              MED15_plus_3 = RNA_ct_yes_MED15_3 / DNA_ct_yes_MED15_3,
                                                              MED15_plus_4 = RNA_ct_yes_MED15_4 / DNA_ct_yes_MED15_4,
                                                              TAF13_minus_1 = RNA_ct_no_TAF13_1 / DNA_ct_no_TAF13_1,
                                                              TAF13_minus_2 = RNA_ct_no_TAF13_2 / DNA_ct_no_TAF13_2,
                                                              TAF13_minus_3 = RNA_ct_no_TAF13_3 / DNA_ct_no_TAF13_3,
                                                              TAF13_minus_4 = RNA_ct_no_TAF13_4 / DNA_ct_no_TAF13_4,
                                                              SPT7_minus_1 = RNA_ct_no_SPT7_1 / DNA_ct_no_SPT7_1,
                                                              SPT7_minus_2 = RNA_ct_no_SPT7_2 / DNA_ct_no_SPT7_2,
                                                              SPT7_minus_3 = RNA_ct_no_SPT7_3 / DNA_ct_no_SPT7_3,
                                                              SPT7_minus_4 = RNA_ct_no_SPT7_4 / DNA_ct_no_SPT7_4,
                                                              MED15_minus_1 = RNA_ct_no_MED15_1 / DNA_ct_no_MED15_1,
                                                              MED15_minus_2 = RNA_ct_no_MED15_2 / DNA_ct_no_MED15_2,
                                                              MED15_minus_3 = RNA_ct_no_MED15_3 / DNA_ct_no_MED15_3,
                                                              MED15_minus_4 = RNA_ct_no_MED15_4 / DNA_ct_no_MED15_4)

## calculate log2-fold change from auxin treatment
UAS_count_norm <- UAS_count_norm %>% ungroup() %>% rowwise() %>% 
  mutate(MED15_minus_mean_a = mean(c(MED15_minus_1, MED15_minus_2), na.rm = TRUE),
         MED15_minus_mean_b = mean(c(MED15_minus_3, MED15_minus_4), na.rm = TRUE),
         MED15_plus_mean_a = mean(c(MED15_plus_1, MED15_plus_2), na.rm = TRUE),
         MED15_plus_mean_b = mean(c(MED15_plus_3, MED15_plus_4), na.rm = TRUE),
         SPT7_minus_mean_a = mean(c(SPT7_minus_1, SPT7_minus_2), na.rm = TRUE),
         SPT7_minus_mean_b = mean(c(SPT7_minus_3, SPT7_minus_4), na.rm = TRUE),
         SPT7_plus_mean_a = mean(c(SPT7_plus_1, SPT7_plus_2), na.rm = TRUE),
         SPT7_plus_mean_b = mean(c(SPT7_plus_3, SPT7_plus_4), na.rm = TRUE),
         TAF13_minus_mean_a = mean(c(TAF13_minus_1, TAF13_minus_2), na.rm = TRUE),
         TAF13_minus_mean_b = mean(c(TAF13_minus_3, TAF13_minus_4), na.rm = TRUE),
         TAF13_plus_mean_a = mean(c(TAF13_plus_1, TAF13_plus_2), na.rm = TRUE),
         TAF13_plus_mean_b = mean(c(TAF13_plus_3, TAF13_plus_4), na.rm = TRUE)) %>%
  mutate(MED15_FC_a = log(MED15_plus_mean_a / MED15_minus_mean_a, 2),
         MED15_FC_b = log(MED15_plus_mean_b / MED15_minus_mean_b, 2),
         SPT7_FC_a = log(SPT7_plus_mean_a / SPT7_minus_mean_a, 2),
         SPT7_FC_b = log(SPT7_plus_mean_b / SPT7_minus_mean_b, 2),
         TAF13_FC_a = log(TAF13_plus_mean_a / TAF13_minus_mean_a, 2),
         TAF13_FC_b = log(TAF13_plus_mean_b / TAF13_minus_mean_b, 2))


# bring in previous dataset with information on coactivator class and gene characteristics
endo_data <- read.csv(file = "//home/jschofie/2303_coactivator_data.csv")
colnames(endo_data) <- c("UAS_name", "UAS_common_name", "UAS_TATA", "UAS_class", "UAS_tail_dep", "UAS_endo_expression", "Spt7_endo", "Taf13_endo", "Med15_endo")

UAS_count_endo_joined <- left_join(UAS_count_norm, endo_data)
UAS_count_endo_joined$core_ID <- factor(UAS_count_endo_joined$core_ID,
                                        levels = c('KRS1','SIT1', 'NOP13', 'RPC10'),ordered = TRUE)

## Pivot table to perform analyses by core promoter
UAS_count_master <- UAS_count_endo_joined %>% ungroup() %>%
  pivot_wider(id_cols = c(UAS_name, core_ID), names_from = core_ID, values_from = c(MED15_minus_mean_a, MED15_plus_mean_a, SPT7_minus_mean_a, SPT7_plus_mean_a,
                                                                                    TAF13_minus_mean_a, TAF13_plus_mean_a, MED15_FC_a, SPT7_FC_a, TAF13_FC_a,
                                                                                    MED15_minus_mean_b, MED15_plus_mean_b, SPT7_minus_mean_b, SPT7_plus_mean_b,
                                                                                    TAF13_minus_mean_b, TAF13_plus_mean_b, MED15_FC_b, SPT7_FC_b, TAF13_FC_b, Spt3_endo, Spt7_endo, Taf1_endo, Taf13_endo, Med15_endo,
                                                                                    UAS_TATA, UAS_class, UAS_tail_dep, UAS_endo_expression))


## Bring in information about control genomic regions
control_assignments <- read.csv(file = "UAS_control_region_assignments.csv")
colnames(control_assignments) <- c("chromosome", "start", "end", "UAS_name", "FAIRE_score", "strand", "type", "join_ID")
control_assignments$UAS_name <- as.character(control_assignments$UAS_name)
UAS_count_master <- left_join(UAS_count_master, control_assignments)
UAS_count_master <- UAS_count_master %>%
  mutate(chromatin = ifelse(type %in% c("AC", "AN"), "accessible",
                            ifelse(type %in% c("IC", "IN"), "inaccessible", "UAS")))


## Bringing in 135bp sequence window surrounding called TFIIB location. Define weak TATA as containing TWTWWA motif upstream of/overlapping center TFIIB location
TFIIB_locs <- read.table(file = "TFIIB_window.txt")
TFIIB_locs <- TFIIB_locs %>% mutate(upstream = substring(V2, 1, 70)) %>%
  mutate(TWTWWA = ifelse(grepl("TTTTTA", upstream) == TRUE, TRUE,
                         ifelse(grepl("TTTTAA", upstream) == TRUE, TRUE,
                                ifelse(grepl("TTTAAA", upstream) == TRUE, TRUE,
                                       ifelse(grepl("TATAAA", upstream) == TRUE, TRUE,
                                              ifelse(grepl("TATTTA", upstream) == TRUE, TRUE,
                                                     ifelse(grepl("TATATA", upstream) == TRUE, TRUE,
                                                            ifelse(grepl("TTTATA", upstream) == TRUE, TRUE,
                                                                   ifelse(grepl("TATTAA", upstream) == TRUE, TRUE, FALSE)))))))))
UAS_count_master <- left_join(UAS_count_master, TFIIB_locs)
UAS_count_master <- UAS_count_master %>%
  mutate(TATA_stat = ifelse(UAS_TATA_KRS1 == "TATA-containing", "strong",
                            ifelse(UAS_TATA_KRS1 == "TATA-less" & TWTWWA == TRUE, "weak", "none")))

### Adding in genes with mediator ChEC-seq signal in UAS (Warfield et al. 2022)
mediator_signal <- read.csv("Mediator_genes_Med8_Med17_signals.csv")
colnames(mediator_signal) <- c("UAS_name", "gene_name", "Med8_sig", "Med17_sig")
UAS_count_master <- left_join(UAS_count_master, mediator_signal)
UAS_count_master <- UAS_count_master %>%
  mutate(Med8_chec = ifelse(!is.na(Med8_sig), TRUE, FALSE),
         Med17_chec = ifelse(!is.na(Med17_sig), TRUE, FALSE)) %>%
  mutate(scaled_express_SIT1 = MED15_minus_mean_a_SIT1*10000,
         scaled_express_KRS1 = MED15_minus_mean_a_KRS1*10000,
         scaled_express_NOP13 = MED15_minus_mean_a_NOP13*10000,
         scaled_express_RPC10 = MED15_minus_mean_a_RPC10*10000)


## Calculate variance across core promoters versus variance across batches, and ratio of activation for CR vs. TFIID core promoters
UAS_count_master <- UAS_count_master %>%
  mutate(ex_mean_a = (MED15_minus_mean_a_KRS1 + MED15_minus_mean_a_SIT1 + MED15_minus_mean_a_NOP13 + MED15_minus_mean_a_RPC10)/4,
         ex_mean_b = (MED15_minus_mean_b_KRS1 + MED15_minus_mean_b_SIT1 + MED15_minus_mean_b_NOP13 + MED15_minus_mean_b_RPC10)/4,
         variance_a = (abs(MED15_minus_mean_a_KRS1 - ex_mean_a) + abs(MED15_minus_mean_a_SIT1 - ex_mean_a) + abs(MED15_minus_mean_a_NOP13 - ex_mean_a) + abs(MED15_minus_mean_a_RPC10 - ex_mean_a))/ex_mean_a,
         variance_b = (abs(MED15_minus_mean_b_KRS1 - ex_mean_b) + abs(MED15_minus_mean_b_SIT1 - ex_mean_b) + abs(MED15_minus_mean_b_NOP13 - ex_mean_b) + abs(MED15_minus_mean_b_RPC10 - ex_mean_b))/ex_mean_b,
         interrep_variance = (abs(MED15_minus_mean_a_KRS1 - MED15_minus_mean_b_KRS1) + abs(MED15_minus_mean_a_SIT1 - MED15_minus_mean_b_SIT1) +
                                abs(MED15_minus_mean_a_NOP13 - MED15_minus_mean_b_NOP13) + abs(MED15_minus_mean_a_RPC10 - MED15_minus_mean_b_RPC10))/((ex_mean_a + ex_mean_b)/2),
         mean_variance = (variance_a + variance_b)/2,
         variance_ratio = mean_variance/interrep_variance,
         variance_class = ifelse(log(variance_ratio, 2) > 1.4, "High",
                                 ifelse(log(variance_ratio, 2) < -1.64, "Low", "Intermediate")),
         class_ratio = log((MED15_minus_mean_a_KRS1 + MED15_minus_mean_b_KRS1 + MED15_minus_mean_b_SIT1 + MED15_minus_mean_a_SIT1)/(MED15_minus_mean_a_NOP13 + MED15_minus_mean_b_NOP13 + MED15_minus_mean_b_RPC10 + MED15_minus_mean_a_RPC10)))


saveRDS(UAS_count_master, paste0("UAS_count_master_", Date, ".rds"))

rm(list = ls())
