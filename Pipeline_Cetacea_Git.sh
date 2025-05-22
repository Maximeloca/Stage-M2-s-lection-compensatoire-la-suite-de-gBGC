####ne garder que les espèces que l'on connait
cd ~/Cetacea/Aligned_Sequences/Align_AA
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align.txt
mkdir align_only_know
cat ~/Cetacea/Aligned_Sequences/Align_AA/list_exons_align.txt | parallel -j50 python3 ~/scripts/script_only_cetace_know.py ~/Cetacea/Aligned_Sequences/Align_AA/{} ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/ ~/Cetacea/list_id_corrected.txt

###remplacer les id par les noms
cd align_only_know
mkdir align_id_name
cat ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/list_exons_align.txt | parallel -j50 python3 ~/scripts/script_changer_nom.py ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/{} ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name/ ~/Cetacea/tab_species_id_corrected.csv


###retirer les espèces qui sont restés en doublon, les dauphins
cd ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align.txt
mkdir align_wt_double
cat ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name/list_exons_align.txt | parallel -j50 python3 ~/scripts/keep-longest_V2.py ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name/{} ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name/align_wt_double/

##
### retirer les séquences avec des codons stop
cd ~/Cetacea/Aligned_Sequences/Align_NT/align_only_know/align_id_name/align_wt_double
find -maxdepth 1 -type f -name "*.fasta" | sed 's/_seq_NT_aligned.fasta//' | cut -f2 -d"/" > list_exons_rm_stop_codon.txt


mkdir wt_stop
cd wt_stop

cat ~/Cetacea/Aligned_Sequences/Align_NT/align_only_know/align_id_name/align_wt_double/list_exons_rm_stop_codon.txt | parallel -j40 python3 ~/scripts/remove_seq_with_stop_codon_cetace.py {} ~/Cetacea/Aligned_Sequences/Align_AA/align_only_know/align_id_name/align_wt_double/{}_seq_AA_aligned.fasta ~/Cetacea/Aligned_Sequences/Align_NT/align_only_know/align_id_name/align_wt_double/{}_seq_NT_aligned.fasta ~/Cetacea/Aligned_Sequences/Align_NT/align_only_know/align_id_name/align_wt_double/wt_stop/ {}_seq_with_stop.csv


mkdir not_used

for align in *.fasta ; do 
if (($(grep -c ">" $align) < 4)) ; then mv $align not_used ; fi 
done

find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align.txt




### PHYLOGENY ################################################################################

## create a file of concatenation of all exons
cd ~/Cetacea

# concatenate exons sequences per specie (one fasta file per species that contains all exons sequences)
mkdir Aligned_concatenated_sequences
cd Aligned_concatenated_sequences

cat ~/Cetacea/list_species.txt | parallel -j50 ~/scripts/script_concatenate_exons_cetace.sh {}

# reformate fasta files
for sp in $(cat ~/Cetacea/list_species.txt); do sed -i ':a;N;/>/!s/\n//;ta;P;D' ${sp}_all_seq_NT_aligned_final.fasta ; done
###permet de suprrimer les retours à la ligne dans les séquences : 


for sp in $(cat ~/Cetacea/list_species.txt); do sed -i 's/!/N/g' ${sp}_all_seq_NT_aligned_final.fasta ; done
###remplace juste les ! par des N dans les séquences 



mkdir wt_chevrons
cd wt_chevrons

for sp in $(cat ~/Cetacea/list_species.txt); do grep -v ">" ../${sp}_all_seq_NT_aligned_final.fasta > ${sp}_all_seq_NT_aligned_final.fasta  ; done

# concatenate sequences of all species in one fasta file
cd ../
mkdir All_seq
cd All_seq

~/scripts/script_concatenate_all_exons_cetace.sh Cetacea

## make the phylogeny on the concatenation of all exons
cd ~/Cetacea
mkdir Phylogeny
cd Phylogeny


# you need to install IQTREE-2


###PHYLOGÉNIE GLOBALE 
iqtree2 -s ~/Cetacea/Aligned_concatenated_sequences/All_seq/Cetacea_all_seq_NT_aligned_final.fasta  -nt 50 -m GTR -bb 1000 2>error_iqtree.txt >out_iqtree.txt
## -s veut dire qu'on va lui donner une séquence 
## -nt les CPU --> mettre que 20 peut ê car BIGMEN bien rempli
## -m --> le modèle 
## GTR --> modèle phylogénétique qu'on a utilisé 
## b --> les bootstraps on a mis 1000

mv ~/Cetacea/Aligned_concatenated_sequences/All_seq/Cetacea_all_seq_NT_aligned_final.fasta.treefile ~/Cetacea/Phylogeny/
### juste iqtree crée l'arbre dans le fichier ou on a la séquence donc on le déplace dans le dossier phylogénie 

cd ~/Cetacea/Phylogeny

python3 ~/scripts/longueur_tot_arbre.py Cetace_multiple_rooted_tree.treefile longueur_arbre.txt

## Root the phylogeny
Rscript ~/scripts/script_multiple_root_tree.R ~/Cetacea/Phylogeny/Cetacea_all_seq_NT_aligned_final.fasta.treefile Sus_scrofa Equus_caballus Vicugna_pacos Camelus_bactrianus Hippopotamus_amphibius Choeropsis_liberiensis Panthalops_hodgsonii Ovis_aries Oryx_leucoryx Gazella_arabica Tragelaphus_eurycerus Bos_taurus Cetace_multiple_rooted_tree.treefile
### on doit voir quelle espèce il faut mettre pour raciner 
### il faut peut être même modifier la ligne de commande pour mettre deux espèces pour raciner car les deux plus eloignées sont des espèces soeurs avec la même distance 
### JE TESTE D'ABORD AVEC UNE SEULE ESPÈCE EN RACINE ET SI JE N'AI PAS LA MÊME TOPOLOGIE IL FAUT CREUSER, JE PEUX REGARDER LES PHYLOS SUR SEAVIEW 
### vérifier la phylo sur l'article d'où proviennent les données 

## make one phylogeny per gene by pruning the global phylogeny
	# keep only the common species between all exons of a same gene
	# on ne garde que les espèces qui sont communes à tous les exons appartenant au même gène
#cd ~/Murinae/Rattini/Aligned_Sequences/alignments/alignments_wt_ref/alignments_wt_ref_filtered/alignments_wt_ref_filtered_wt_stop/

#find -maxdepth 1 -type f -name "*.fasta" | cut -f1,2,3 -d"_" | cut -f2 -d"/" > Genes_exons_pos_list.csv
##on crée un fichier csv qui contient gene_exon_position
#sed -i 's/_/\t/g' Genes_exons_pos_list.csv
## remplace le séparateur qui est _ par un tab

#cut -f1 Genes_exons_pos_list.csv | sort -u > Genes_list.txt ## liste de gènes
cd ~/Cetacea/Aligned_Sequences
mkdir FINAL
cd FINAL

for file in ~/Cetacea/Aligned_Sequences/Align_NT/align_only_know/align_id_name/align_wt_double/wt_stop/*.fasta ; do
    cp "$file" ~/Cetacea/Aligned_Sequences/FINAL/
done

### ne garder que les espèces avec plus de 4 espèces
mkdir not_used
for align in *.fasta_seq_NT_aligned.fasta ; do 
if (($(grep -c ">" $align) < 4)) ; then mv $align not_used ; fi 
done
###on garde que les séquences comptant au minimum 4 espèces 

####On ne garde que les séquences avec une longueur de plus de 50 sites

##test
#cd ~/Cetacea/Aligned_Sequences/30
#find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons.txt
#cat ~/Cetacea/Aligned_Sequences/30/list_exons.txt | parallel -j50 python3 ~/scripts/leng_align.py ~/Cetacea/Aligned_Sequences/30/{} ~/Cetacea/Aligned_Sequences/30/

#cat ~/Cetacea/Aligned_Sequences/30/list_exons.txt | parallel -j50 ~/scripts/mv_leng_inf_50.sh {} ~/Cetacea/Aligned_Sequences/30/ ~/Cetacea/Aligned_Sequences/30/inf_50

####On ne garde que les séquences avec une longueur de plus de 50 sites
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align_sup_4.txt

#on crée une liste des fichiers avec mla longueur des séquences
cat ~/Cetacea/Aligned_Sequences/FINAL/list_exons_align_sup_4.txt | parallel -j50 python3 ~/scripts/leng_align.py ~/Cetacea/Aligned_Sequences/FINAL/{} ~/Cetacea/Aligned_Sequences/FINAL

##On déplace les fichiers ayant moins de 50 bp
mkdir inf_50
cat ~/Cetacea/Aligned_Sequences/FINAL/list_exons_align_sup_4.txt | parallel -j50 ~/scripts/mv_leng_inf_50.sh {} ~/Cetacea/Aligned_Sequences/FINAL/ ~/Cetacea/Aligned_Sequences/FINAL/inf_50


####on retire les outgroup
#test
#find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align_sup_4_sup_50.txt
#cat ~/Cetacea/Aligned_Sequences/30/list_exons_align_sup_4_sup_50.txt | parallel -j50 python3 ~/scripts/remove_out_group.py ~/Cetacea/Aligned_Sequences/30/{} ~/Cetacea/list_outgoup_test.txt ~/Cetacea/Aligned_Sequences/30/wt_outgroup/

cd ~/Cetacea/Aligned_Sequences/FINAL 
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_align_sup_4_sup_50.txt
mkdir wt_outgroup
cat ~/Cetacea/Aligned_Sequences/FINAL/list_exons_align_sup_4_sup_50.txt | parallel -j50 python3 ~/scripts/remove_out_group.py ~/Cetacea/Aligned_Sequences/FINAL/{} ~/Cetacea/list_outgroup.txt ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/


cd ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_final.txt

###on retire les séquences avec moins de 4 espèces car maintenant on a enlevé les outgroup
mkdir not_used
for align in *.fasta_seq_NT_aligned.fasta ; do 
if (($(grep -c ">" $align) < 4)) ; then mv $align not_used ; fi 
done

cd ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup
find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_final.txt


cd ~/Cetacea/Phylogeny/
mkdir exon_tree
cd exon_tree

###retirer aléatoirement dans les alignements ayant plus de 25 esp, des sp jusqu'à arriver à 25

#compter le nb d'espèces par alignements 
cat ~/Cetacea2/Aligned_Sequences/FINAL/list_exons_final.txt | parallel -j100 ~/scripts/count_sp_align.sh ~/Cetacea2/Aligned_Sequences/FINAL/{}


## création d'un arbre par exon
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j50 Rscript ~/scripts/script_prune_tree.R ~/Cetacea/Phylogeny/Cetace_multiple_rooted_tree.treefile ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/{} {}_tree.treefile 2> error_gene_tree.txt >out_gene_tree.txt

### SUPPLEMENTARY FILTER OF PARALOGS #########################################################
### les duplications peuvent modifier longueur de branches 

cd ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup
for fasta in ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/*.fasta; do sed -i 's/!/N/g' "$fasta" ; done
###remplace juste les ! par des N dans les séquences 

for file in *.fasta ; do len=$(seqkit fx2tab $file -n -l | head -1 | cut -f2) ; echo -e "${file}\t${len}" ; done > align_len.csv

cd ~/Cetacea/
mkdir Sup_filter_paralogs
cd Sup_filter_paralogs

   ## !! ## chek that you have the python module ete3 installed ##
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j50 ~/scripts/script_filtre_sup_paralogues_cetaces.sh {} > tab_lg_tree_for_sup_filter_paralogues.csv ###attention vu qu'on a pas le gène on a enlevé dans le tableau tab_lg_tree_for_sup_filter_paralogues la colonne gèhe

#savoir si les fichiers d'erreurs sont vides ou non --> ici 
find -name "error_*.txt" -exec sh -c 'if [ -s "$1" ]; then echo "$1 is not empty"; else echo "$1 is empty"; fi' _ {} \; > output_log.txt


~/scripts/script_list_paralog_exons_cetace.R tab_lg_tree_for_sup_filter_paralogues.csv ###attention vu qu'on a pas le gène on a enlevé dans le tableau tab_lg_tree_for_sup_filter_paralogues la colonne gèhe, il faut donc modifier le script R aussi 
### crée le fichier paralog_exons.txt

###vérifier combien d'exon il nous reste après retirer paralogue 



### SUBSTITUTIONS MAPPING ####################################################################


cd ~/Cetacea/
mkdir Substitutions_mapping





#cut -f2 Genes_exons_nbsp_min10.csv >  list_exons_min10sp.txt

### make the substitutions mapping

cd ~/Cetacea/Substitutions_mapping/
mkdir output

# you need to install bppml and mapnh (put them in your bin)

#cat ~/Murinae/Rattini/Aligned_Sequences/FINAL/list_exons_min10sp.txt | parallel -j50 ~/scripts/script_mapNH_bpp.sh Rattini {} 2> error_sortie_mapNH.txt > out_mapNH.txt


###garder list avec au dessus de 5 sp
###bppml permet de mesurer les longeurs de branches de l'arbres
####à la fin de cette commande on aura les 8 fichiers avec les substitutions

###même commande que juste au dessus mais en utilisant pas les exons ayant plus de 10 espèces donc avec la liste Final_exons_list.txt
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j50 ~/scripts/script_mapNH_cetace_bpp.sh Cetacea {} 2> error_sortie_mapNH.txt > out_mapNH.txt


#~/scripts/script_mapNH_bpp.sh Rattini ENSMUSE00000284194 2> error_sortie_mapNH.txt > out_mapNH.txt


find -type f -name "*.dnd_1" | cut -f1 -d"_" | cut -f2 -d"/" > list_exons_mapping_ok.txt
find -type f -name "*.dnd_1" | sed 's/_ml\.dnd_1$//' | cut -f2 -d"/" > list_exons_mapping_ok.txt



cd output
mkdir correct_file_names
cat ../list_exons_mapping_ok.txt | parallel -j50 ~/scripts/transform_mapping_file_name.sh {} ### on change juste le nom de fichier et on le transforme en .csv


####sortir le plot des substitutions
cd correct_file_names
mkdir distrib_subst
cd distrib_subst
#cat ~/Murinae/Aligned_Sequences/Hydromyini/FINAL/list_exons_min10sp.txt | parallel -j20 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Murinae/Substitutions_mapping/Hydromyini/output/correct_file_names/
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j40 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Cetacea/Substitutions_mapping/output/correct_file_names/

#on concatène tous les fichiers 
for file in *_in_branches.csv ; do cat $file ; done > count_subst_types_in_branches_all_exons.csv

### DETECTION OF gBGC EPISODES ####################################################################

#######
## 1 ## Calculate statistics ###
#######

### Distribution of synonymous substitutions

cd ~/Cetacea
mkdir Detecting_gBGC
cd Detecting_gBGC
mkdir Distrib_subst_S
cd Distrib_subst_S

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j50 Rscript ~/scripts/script_nb_subst.R {} ~/Cetacea/Substitutions_mapping/output/correct_file_names/ {}_stat_nb_subst.csv 2> errors_script_nb_subst.txt > out_script_nb_subst.txt

#somme des substitutions pour tous les codons d'une branche


#######
## 2 ## Randomisations of statistics ###
#######

cd ~/Cetacea
mkdir Randomisations
cd Randomisations
mkdir Br_lg_tables
mkdir Distrib_subst
mkdir I_Moran

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j50 ~/scripts/script_simul_stat_cetace.sh {} Cetacea 1000 > out_simul.txt 2> error_simul.txt


mkdir used
for exon in $(cat ~/Murinae/Aligned_Sequences/Rattini/FINAL/list_exons_sup1000bp.txt) ; do #cat ~/Murinae/Rattini/Aligned_Sequences/FINAL/list_exons_sup1000bp.txt
mv ${exon}_I_simul.csv used ;
done
rm *.csv
cd used
mv *.csv ..
cd ..
rm -r used

#######
## 3 ## Detection of gBGC episodes ###
#######

cd ~/Cetacea//Detecting_gBGC/
mkdir Detection
cd Detection
mkdir signif_stat
cd signif_stat

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j20 ~/scripts/script_signif_stat_cetace.sh {} Cetacea 1000 >out_commande_signif_stat.txt 2> error_commande_signif_stat.txt

#for ex in $(cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt) ; do len=$(grep ${ex} ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/exons_len.csv | cut -f2) ; gene=$(grep ${ex} ~/Cetacea/Reference_exons/Gene_exons_ref_list.txt | cut -f1) ; pos=$(grep -w ${ex} ~/Murinae/Reference_exons/Gene_exons_ref_list.txt | cut -f3) ; echo -e "${gene}\t${ex}\t${pos}\t${len}" ; done > gene_exon_lenseq_mapping_ok.txt
for ex in $(cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt) ; do len=$(grep ${ex} ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/exons_len.csv | cut -f2) ; echo -e "${ex}\t${len}" ; done > exon_lenseq_mapping_ok.txt


cd ~/Cetacea/Detecting_gBGC/Detection/
mkdir episodes
cd episodes

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j50 ~/scripts/script_detection_episodes_cetace.sh {} Cetacea 1000 995 > out_detection_episodes.txt 2> errors_detection_episodes.txt



# Parcours des fichiers *_episodes_detection.csv
for file in *_episodes_detection.csv; do
    # Comptage des occurrences
    nb_yes=$(grep -c -w "YES" "$file")
    nb_sureno=$(grep -c -w "SURENO" "$file")
    nb_no=$(grep -c -w "NO" "$file")

    # Affichage des résultats dans un fichier
    echo -e "${file}\t${nb_yes}\t${nb_sureno}\t${nb_no}" >> tab_episode.csv
done
#    count_no="$count_no" + "$nb_no"
#    count_sureno="$count_sureno" + "$nb_sureno"
#    count_yes="$count_yes" + "$nb_yes"

### DISTRIBUTION & RANDOMISATION OF NON-SYNONYMOUS SUBSTITUTIONS ####################################

#######
## 1 ## Distribution ###
#######

cd ~/Cetacea/
mkdir Impact_episodes_NS
cd Impact_episodes_NS

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j10 Rscript ~/scripts/script_nb_subst_NS.R {} ~/Cetacea/Substitutions_mapping/output/correct_file_names/ {}_stat_nb_subst_NS.csv 2> errors_script_nb_subst_NS.txt > out_script_nb_subst_NS.txt

#######
## 2 ## Randomisation ###
#######

cd ~/Cetacea/Randomisations/
mkdir Distrib_NS
cd Distrib_NS

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j40 ~/scripts/script_simul_distrib_subst_NS_cetace.sh {} Cetacea 1000 > sortie_simul_NS.txt 2> error_simul_NS.txt

### CALCULATE GC CONTENT #########################################################################


cd ~/Cetacea
mkdir GC_content
cd GC_content

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j10 ~/scripts/script_GC_content_cetace.sh {} Cetacea >out_GC_content.txt 2>errors_GC_content.txt

### BRANCH RELATIONS #############################################################################
####il faut savoir comment les branches sont numérotés dans les arbres et ici c'est du postorder

#cd ~/Cetacea/Substitutions_mapping/

#for exon in $(cat list_exons_mapping_ok.txt) ; do
#gene=$(grep $exon ~/Murinae/Reference_exons/Gene_exons_ref_list.txt | cut -f1) ; echo -e "${gene}\t${exon}" ; done > list_genes_exons_mapping_ok.txt

#cut -f1 list_genes_exons_mapping_ok.txt | sort -u > list_genes_mapping_ok.txt

cd ~/Cetacea/Detecting_gBGC/Detection/
mkdir relations_branches
cd relations_branches

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j50 ~/scripts/script_num_nodes_tree_cetace.sh {} Cetacea 2> error_rel_br.txt > out_rel_br.txt
#crée un tableau avec 4 colonne pour chaque gène: 
# numéro_de_branche		numéro_branche_ascendante 	numéro_branche_descendante_1	numéro_branche_descendante_2

### FINAL TABLE WITH INFORMATIONS FOR ALL EXONS FOR GENERAL RESULTS ###############################

cd ~/Cetacea/
mkdir Final_data
cd Final_data

cat ~/Cetacea/Substitutions_mapping/list_exons_mapping_ok.txt | parallel -j50 python3 ~/scripts/script_table_recap_cetace.py {} ~/Cetacea/Substitutions_mapping/exon_lenseq_mapping_ok.txt ~/Cetacea/Detecting_gBGC/Detection/relations_branches/{}_rel_br.csv ~/Cetacea/Randomisations/Br_lg_tables/{}_br_len.csv ~/Cetacea/Detecting_gBGC/Detection/episodes/{}_episodes_detection.csv ~/Cetacea/Detecting_gBGC/Distrib_subst_S/{} ~/Cetacea/Detecting_gBGC/Detection/signif_stat/{} ~/Cetacea/Impact_episodes_NS/{} ~/Cetacea/GC_content/{} {}_tab_recap.csv > sortie_tab_recap.txt 2> error_tab_recap.txt

for file in *_recap.csv ; do sed 1d $file ; done > tab_recap_all_genes.csv

#ls *_recap.csv | cut -f1 -d'_' > list_genes.txt

### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Cetacea/Final_data
find . -maxdepth 1 -type f -name "*_tab_recap.csv" | sed 's/_tab_recap.csv//' | cut -f2 -d"/" > list_exons_final.txt

cd ~/Cetacea/
mkdir Final_data_analyse_NonSynonymous_Compensation
cd Final_data_analyse_NonSynonymous_Compensation

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j50 python3 ~/scripts/analyse_random_NS_cetace.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh

####CRÉATION D'UN FILTRE PERMETTANT DE VOIR LA QUALITÉ DES ALIGNEMENTS 
cd ~/Cetacea/
mkdir filtre_align
cd filtre_align

python3 ~/scripts/script_list_exons_algnt_error.py ~/Cetacea/Final_data/list_exons_final.txt ~/Cetacea/Final_data/ 0.05 2>error.txt >out.txt

###A partir de 0,1 on dirait que ca devient bizarre je crois, continuer à regarder ceux la

cat ~/Cetacea/filtre_align/list_exons_0.01 | parallel -j40 ~/scripts/cp_fichier_cetace.sh {} 0.01

###### SORTING EPISODE TYPE (punctual, two-branches, or all-exons) ####################################

cd ~/Cetacea/
mkdir Sort_episodes_types
cd Sort_episodes_types

### make table for the sorting ###
#cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j30 python3 ~/scripts/script_tab_for_sort_ep_cetace.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt 2>error.txt >out.txt
###ce script récupère dans les gènes possédant au moins un exon avec un épisode des informations sur la branche avec l'épisode pour tous les exons donc même ceux qui n'ont pas d'épisodes et les branches descendantes à la branche contenant l'épisode 
###un script certainement pour voir faire le tri dans les épisodes et les trier 


###crée le tableau des exons et des cds associées
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j40 ~/scripts/tab_exon_cds.sh {}


cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j50 python3 ~/scripts/script_tab_for_sort_ep_cetace_V2.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt ~/Cetacea/Sort_episodes_types_wt_gene/tab_exon_cds_header.csv 2>error.txt >out.txt


for file in *_tab_for_sort_ep.csv ; do sed 1d $file ; done > tab_for_sort_ep_all_genes_Cetacea_wt_header.csv

## AJOUTER TITRE COLONNES ###
  # ajout du header au tableau
echo ""Gene,Exon,GC3,Lg_seq,Br_ep,Episode,Nb_S_WS_ep,Nb_S_SW_ep,Nb_S_WS_postep,Nb_S_SW_postep"" > tab_for_sort_ep_all_genes_Cetacea.csv
cat tab_for_sort_ep_all_genes_Cetacea_wt_header.csv >> tab_for_sort_ep_all_genes_Cetacea.csv


find -maxdepth 1 -type f -name "*_tab_for_sort_ep.csv" | cut -f2 -d"/" | sed 's/_tab_for_sort_ep.csv$//' > list_ep.csv



### make the sorting ####
~/bin/sort_ep_type/MLsort_single_exons tab_for_sort_ep_all_genes_Cetacea.csv optfile.txt posterior_optimum_Cetacea.csv > output




####relancer cette partie en ne mettant qu'un seul fichier de sortie pour identifier les gènes avec les erreurs 
mkdir Data2
cd Data2
cat ~/Murinae/Rattini/Final_data/list_genes.txt | parallel -j40 python3 ~/scripts/dist_comp.py {} ~/Murinae/Rattini/Final_data/ ~/Murinae/Rattini/Sup_filter_paralogs/paralog_exons.txt ~/Murinae/Rattini/Substitutions_mapping/output/correct_file_names/ &> out_dist_comp_error.txt
python3 ~/scripts/dist_comp_test.py ENSMUSG00000078768 ~/Murinae/Rattini/Final_data/ ~/Murinae/Rattini/Sup_filter_paralogs/paralog_exons.txt ~/Murinae/Rattini/Substitutions_mapping/output/correct_file_names/ &> out_dist_comp_error.txt
python3 ~/scripts/dist_comp_test.py ENSMUSG00000000167 ~/Murinae/Rattini/Final_data/ ~/Murinae/Rattini/Sup_filter_paralogs/paralog_exons.txt ~/Murinae/Rattini/Substitutions_mapping/output/correct_file_names/ &> out_dist_comp_error.txt
python3 ~/scripts/dist_comp_test.py ENSMUSG00000018168 ~/Murinae/Rattini/Final_data/ ~/Murinae/Rattini/Sup_filter_paralogs/paralog_exons.txt ~/Murinae/Rattini/Substitutions_mapping/output/correct_file_names/ &> out_dist_comp_error.txt





 ####CRÉATION D'UN FILTRE PERMETTANT DE VOIR LA QUALITÉ DES ALIGNEMENTS 
 cd ~/Murinae/Rattini
 mkdir filtre_align
 cd filtre_align

python3 ~/scripts/script_list_exons_algnt_error.py ~/Murinae/Rattini/Final_data/list_genes.txt ~/Murinae/Rattini/Final_data/ 0.10 2>error.txt >out.txt

###A partir de 0,1 on dirait que ca devient étrange je crois, continuer à regarder ceux la

cat ~/Murinae/Rattini/filtre_align/list_exons_0.15 | parallel -j40 ~/scripts/cp_fichier.sh {} 0.15

cd ~/Cetacea/Substitutions_mapping
mkdir test_wt_extrem_gap
cd test_wt_extrem_gap
mkdir test_wt_extrem_gap align_wt_gap

cd align_wt_gap 
cat ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/align_wt_gap/list_exons_final.txt | parallel -j70 ~/scripts/script_mapNH_cetace_bpp_test.sh Cetacea {} align_wt_gap 2> error_sortie_mapNH.txt > out_mapNH.txt




find -type f -name "*.dnd_1" | sed 's/_ml\.dnd_1$//' | cut -f2 -d"/" > list_exons_mapping_ok.txt



cd output
mkdir correct_file_names
cat ../list_exons_mapping_ok.txt | parallel -j50 ~/scripts/transform_mapping_file_name.sh {} ### on change juste le nom de fichier et on le transforme en .csv


####sortir le plot des substitutions
cd correct_file_names
mkdir distrib_subst
cd distrib_subst
#cat ~/Murinae/Aligned_Sequences/Hydromyini/FINAL/list_exons_min10sp.txt | parallel -j20 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Murinae/Substitutions_mapping/Hydromyini/output/correct_file_names/
cat ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/align_gap/list_exons_final.txt | parallel -j40 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/extreme_gap/output/correct_file_names/
cat ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/align_wt_gap/list_exons_final.txt | parallel -j40 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/wt_gap/output/correct_file_names/

#on concatène tous les fichiers 
for file in *_in_branches.csv ; do cat $file ; done > count_subst_types_in_branches_all_exons_sans_gap.csv



cd ~/Cetacea/Substitutions_mapping
mkdir test_ungap_0.2
cd test_ungap_0.2
cd ungap

find -maxdepth 1 -type f -name "*.fasta" | cut -f2 -d"/" > list_exons_final.txt


cd ..
cat ~/Cetacea/Substitutions_mapping/test_ungap_0.2/ungap/list_exons_final.txt | parallel -j70 ~/scripts/script_mapNH_cetace_bpp_test.sh Cetacea {} ungap 2> error_sortie_mapNH.txt > out_mapNH.txt


cat ~/Cetacea/Substitutions_mapping/test_ungap_0.2/ungap/list_exons_final.txt | parallel -j40 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Cetacea/Substitutions_mapping/test_wt_extrem_gap/extreme_gap/output/correct_file_names/

find -type f -name "*.dnd_1" | sed 's/_ml\.dnd_1$//' | cut -f2 -d"/" > list_exons_mapping_ok.txt



cd output
mkdir correct_file_names
cat ../list_exons_mapping_ok.txt | parallel -j50 ~/scripts/transform_mapping_file_name.sh {} ### on change juste le nom de fichier et on le transforme en .csv


####sortir le plot des substitutions
cd correct_file_names
mkdir distrib_subst
cd distrib_subst
#cat ~/Murinae/Aligned_Sequences/Hydromyini/FINAL/list_exons_min10sp.txt | parallel -j20 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Murinae/Substitutions_mapping/Hydromyini/output/correct_file_names/
cat ~/Cetacea/Substitutions_mapping/test_ungap_0.2/ungap/list_exons_final.txt | parallel -j40 ~/scripts/script_distrib_ratio_WS_SW_subst.R {} ~/Cetacea/Substitutions_mapping/test_ungap_0.2/output/correct_file_names/

#on concatène tous les fichiers 
for file in *_in_branches.csv ; do cat $file ; done > count_subst_types_in_branches_all_exons_sans_gap_0.2.csv



nb_esp_tot=0
for align in *.fasta_seq_NT_aligned.fasta ; do 
    nb_esp=(grep )
if (($(grep -c ">" $align) < 4)) ; then mv $align not_used ; fi 
done








####Partie ou on met les ep NS à 0.8 et les no NS à 0.2

### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Cetacea/
mkdir Final_data_analyse_NonSynonymous_Compensation_2
cd Final_data_analyse_NonSynonymous_Compensation_2

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/analyse_random_NS_cetace_V2.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh




###### SORTING EPISODE TYPE (punctual, two-branches, or all-exons) ####################################

cd ~/Cetacea/
mkdir Sort_episodes_types_2
cd Sort_episodes_types_2

### make table for the sorting ###
#cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j30 python3 ~/scripts/script_tab_for_sort_ep_cetace.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt 2>error.txt >out.txt
###ce script récupère dans les gènes possédant au moins un exon avec un épisode des informations sur la branche avec l'épisode pour tous les exons donc même ceux qui n'ont pas d'épisodes et les branches descendantes à la branche contenant l'épisode 
###un script certainement pour voir faire le tri dans les épisodes et les trier 


###crée le tableau des exons et des cds associées
cat ~/Cetacea/Aligned_Sequences/FINAL/wt_outgroup/list_exons_final.txt | parallel -j100 ~/scripts/tab_exon_cds.sh {}
cp tab_exon_cds.csv tab_exon_cds_wt_header.csv
rm tab_exon_cds.csv


echo -e "Exon\tGene" > tab_exon_cds.csv
cat tab_exon_cds_wt_header.csv >> tab_exon_cds.csv


cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/script_tab_for_sort_ep_cetace_V2.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt ~/Cetacea/Sort_episodes_types_2/tab_exon_cds.csv 2>error.txt >out.txt




for file in *_tab_for_sort_ep.csv ; do sed 1d $file ; done > tab_for_sort_ep_all_genes_Cetacea_wt_header.csv

## AJOUTER TITRE COLONNES ###
  # ajout du header au tableau
echo ""Gene,Exon,GC3,Lg_seq,Br_ep,Episode,Nb_S_WS_ep,Nb_S_SW_ep,Nb_S_WS_postep,Nb_S_SW_postep"" > tab_for_sort_ep_all_genes_Cetacea.csv
cat tab_for_sort_ep_all_genes_Cetacea_wt_header.csv >> tab_for_sort_ep_all_genes_Cetacea.csv

### make the sorting ####

touch optfile.txt
~/bin/sort_ep_type/MLsort_single_exons tab_for_sort_ep_all_genes_Cetacea.csv optfile.txt posterior_optimum_Cetacea.csv > output



####Partie ou on met les ep NS à 0.6 et les no NS à 0.4

### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Cetacea/
mkdir Final_data_analyse_NonSynonymous_Compensation_3
cd Final_data_analyse_NonSynonymous_Compensation_3

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/analyse_random_NS_cetace_V3.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh


### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Cetacea/
mkdir Final_data_analyse_NonSynonymous_Compensation_4
cd Final_data_analyse_NonSynonymous_Compensation_4

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/analyse_random_NS_V5.py {} ~/Cetacea/Final_data/ ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh


###regarder la distribution des substitutions NS

cd Cetacea
mkdir Count_subst_NS
cd Count_subst_NS

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/distrib_subst_br_post_ep.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

for file in *_tab_NS.csv ; do sed 1d $file ; done > tab_all_ep_Cetacea_wt_header.csv




####Refaire avec le bon script 

### FINAL TABLES FOR COMPENSATION ANALYSES ########################################################

cd ~/Cetacea/
mkdir Final_data_analyse_NonSynonymous_Compensation_new
cd Final_data_analyse_NonSynonymous_Compensation_new

cat ~/Cetacea/Final_data/list_exons_final.txt | parallel -j100 python3 ~/scripts/analyse_random_NS_cetace_Marie.py {} ~/Cetacea/Final_data/{} ~/Cetacea/Randomisations/Distrib_NS/{} ~/Cetacea/Sup_filter_paralogs/paralog_exons.txt >out_analyse_randomNS.txt 2>error_analyse_randomNS.txt

~/scripts/script_concat_tab_analyse_randomNS_V3.sh





####mesurer le Dn/Ds

~/scripts/script_mapNH_cetacea_Dn_Ds.sh Cetacea Cetacea_all_seq_NT_aligned_final.fasta 2> error_sortie_mapNH.txt > out_mapNH.txt
 