## Change dir where files live 
cd /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/00_from_condor/

## Copy over files from chtc to Beast server
rsync -av \
--ignore-existing \
--include '*.tar.gz' \
--exclude '*.out' \
--exclude '*.err' \
--exclude '*.txt' \
--exclude '*.log' \
--exclude '*.out' \
--exclude 'test_script.R' \
--exclude 'my_packages.tar.gz' \
--exclude 'packages.tar.gz' \
--exclude 'R402.tar.gz' \
<chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/0* \
/home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/00_from_condor/

## Create a list of files which were transfered onto Beast server
ls $PWD > list.txt

## Pull that back to chtc
scp -r <chtc-user>@<beast-host>:/home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/00_from_condor/list.txt ./

## Delete those files from chtc
xargs -a list.txt -d'\n' rm


# ## Data 10 (Not going to include)
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/10_Seedling_V1_602genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/10_Seedling_V1_602genos* ./
# ## Copy ALL files from condor just to have a backup
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ## Untar all files
# ls *.gz | xargs -n1 tar -xzvf
# ## Put all the text files into a common directory
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/10_Seedling_V1_602genos/ ; done
# ## Clean up the directory
# rm *.tar.gz
# rm -rf 10_Seedling_V1_602genos_chunk*

# ## Data 9
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/09_Seedling_V1_304genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# # for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/ ; done
# rm *.tar.gz
# rm -rf 09_Seedling_V1_304genos_chunk*

# ## Data 8
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/08_Seedling_V1_187genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/08_Seedling_V1_187genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/08_Seedling_V1_187genos/ ; done
# rm *.tar.gz
# rm -rf 08_Seedling_V1_187genos_chunk*

# ## Data 7
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/07_LMAD_131genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/07_LMAD_131genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/07_LMAD_131genos/ ; done
# rm *.tar.gz
# rm -rf 07_LMAD_131genos_chunk*

# ## Data 6
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/06_Kern_169genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/06_Kern_169genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/06_Kern_169genos/ ; done
# rm *.tar.gz
# rm -rf 06_Kern_169genos_chunk*

# ## Data 5
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/05_Gshoot_178genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/05_Gshoot_178genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/05_Gshoot_178genos/ ; done
# rm *.tar.gz
# rm -rf 05_Gshoot_178genos_chunk*

# ## Data 4
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/04_L3Base_178genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/04_L3Base_178genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/04_L3Base_178genos/ ; done
# rm *.tar.gz
# rm -rf 04_L3Base_178genos_chunk*

# ## Data 3
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/03_L3Tip_180genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/03_L3Tip_180genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/03_L3Tip_180genos/ ; done
# rm *.tar.gz
# rm -rf 03_L3Tip_180genos_chunk*

# ## Data 2
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/02_LMAN_183genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/02_LMAN_183genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/02_LMAN_183genos/ ; done
# rm *.tar.gz
# rm -rf 02_LMAN_183genos_chunk*

# ## Data 1
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/01_GRoot_175genos/
# scp -r <chtc-user>@submit-1.chtc.wisc.edu:/home/<chtc-user>/01_GRoot_175genos* ./
# cp * /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/00_from_condor
# ls *.gz | xargs -n1 tar -xzvf
# for f in **/*.txt; do mv "$f" /home/<chtc-user>/eqtl_pdi_coexpr/results/04_eQTL_results_phenolics_and_tfs/01_GRoot_175genos/ ; done
# rm *.tar.gz
# rm -rf 01_GRoot_175genos_chunk*





