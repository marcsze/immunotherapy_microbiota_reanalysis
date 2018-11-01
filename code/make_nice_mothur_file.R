## Script to create a custom guide file for mothur processing
## Marc Sze

library(tidyverse)

# Read in the needed data
sra_data <- read_tsv("data/process/metadata/sra_metadata.txt") %>% 
  select(Library_Name, Sample_Name, tissue, Run)

mothur_file <- read_tsv("data/process/stability.paired.files", col_names = F) %>% 
  select(-X4)


# Create new stability.files

new_mothur_file <- mothur_file %>% 
  left_join(sra_data, by = c("X1" = "Run")) %>% 
  select(Library_Name, X2, X3)

# Write out the new file
write_tsv(new_mothur_file, "data/process/stability.files")
