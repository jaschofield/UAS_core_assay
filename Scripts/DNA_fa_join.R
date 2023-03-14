### R script to join R1 and R2 from DNA-seq into single string
library(tidyverse)

DNA_R1 <- as.data.frame(readLines('sample_R1_001.fastq'))
DNA_R2 <- as.data.frame(readLines('sample_R2_001.fastq'))

head(DNA_R1)

DF_R1 <- as.data.frame(matrix(DNA_R1[,1], byrow=TRUE, ncol = 4))
DF_R2 <- as.data.frame(matrix(DNA_R2[,1], byrow=TRUE, ncol = 4))

head(DF_R1)

DF_R1 <- DF_R1 %>%
separate(V1, c("tag", "indexing"), sep = ' ')

colnames(DF_R1) <- c('tag', 'indexing', 'R1', 'a', 'b')

head(DF_R1)

DF_R2 <- DF_R2 %>%
separate(V1, c("tag", "indexing"), sep = ' ')

colnames(DF_R2) <- c('tag', 'indexing', 'R2', 'a', 'c')

joined_df <- left_join(DF_R1, DF_R2, by = 'tag') %>% filter(!is.na(R2)) %>% mutate(read = paste0(R1, R2))
head(joined_df)
sequences <- joined_df$read
write.table(sequences, file = "sample_DNA.txt")
