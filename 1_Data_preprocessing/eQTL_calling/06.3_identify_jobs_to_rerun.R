setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/")

files_run = as_tibble(list.files("00_from_condor/")) %>% 
  mutate(i = as.numeric(str_extract(value, "\\d+")),
         value1 = as.numeric(str_extract(value, "\\d+(?=\\.chunk)")),
         value2 = as.numeric(str_extract(value, "\\d+(?=\\.tar.gz)"))) %>% 
  filter(!value == "list.txt")
         

sudo_chrom = 150
expr_split = 40
n = 8

chunk1 = as_tibble(rep(seq(1,expr_split,1),sudo_chrom)) %>% 
  group_by(value) %>% 
  mutate(value1 = rep(seq_len(n())))
chunk_list = cbind(chunk1[rep(1:nrow(chunk1), n),], i = rep(c(9,7,6,5,4,3,2,1), each = nrow(chunk1))) %>% 
  select(i,value1,value) %>% 
  setNames(., c("i","value1","value2"))

my_check = full_join(files_run, chunk_list)
my_re_do = my_check[!complete.cases(my_check),2:4]

write.table(my_re_do, "redo_list.txt", col.names = FALSE, row.names = FALSE)

## Pull that back to chtc
scp -r <chtc-user>@<beast-host>:/home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/redo_list.txt ./