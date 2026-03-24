#!/bin/bash

######################################################
# Parameter Passing
while getopts "p:i:h" opt; do
  case $opt in
    p)
      chrpairs=$OPTARG
      ;;
    i)
      max_iqr=$OPTARG
      ;;
    h)
     echo ""
      ;;
#    \?)
#     echo "Invalid option: -$OPTARG"
#     ;;
  esac
done

######################################################
# Help documentation

 display_usage() {
         echo -e "\nThis script is for filtering allele pairs\n"
         echo -e "Usage:  bash allele_filtering.sh -p [] -i []"
         echo -e "Example: bash allele_filtering.sh -p chrpairs.txt -i 3\n"
         echo -e "-p: Necessary parameter. homologous chromosome pairs"
         echo -e "-i: Necessary parameter. the maximum multiple of the IQR (will run for 1, 1.5, 2, ... up to this value with step 0.5, must be >= 1)"
         echo -e "-h or --help: Usage\n"
         }


 # if less than one arguments supplied, display usage
         if [  $# -le 0 ]
         then
                 display_usage
                 exit 1
         fi

 # check whether user had supplied -h or --help . If yes display usage
         if [[ ( $* == "--help") ||  $* == "-h" ]]
         then
                 display_usage
                 exit 0
         fi

######################################################
# Variable Determination

 if [ ! $chrpairs ]|| [ ! $max_iqr ];then
         echo -e "\nERROR: Missing necessary files."
         exit 1
 fi

# Check if max_iqr is a valid numeric value
if ! [[ "$max_iqr" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        echo -e "\nERROR: -i must be a numeric value.\n"
        exit 1
fi

# Check if max_iqr is >= 1
if (( $(echo "$max_iqr < 1" | bc -l) )); then
        echo -e "\nERROR: -i must be >= 1.\n"
        exit 1
fi

#######################################################
# Work Path
pwd=`pwd`

#######################################################
# Generate IQR value sequence (from 1 to max_iqr, step 0.5)
iqr_values=()
current=1
while (( $(echo "$current <= $max_iqr" | bc -l) )); do
    # Format: remove trailing .0 for integers
    clean_current=$(echo "$current" | sed 's/\.0$//')
    iqr_values+=("$clean_current")
    current=$(echo "$current + 0.5" | bc -l)
done

#######################################################
# Run filtering and plotting for each IQR multiple
for iqr in "${iqr_values[@]}"; do
    echo -e "\n========================================"
    echo -e "Processing IQR multiple: $iqr"
    echo -e "========================================\n"

    # Create separate directory for each multiple
    work_dir="$pwd/02_wgdi-Tukey_${iqr}"
    mkdir -p "$work_dir"
    cd "$work_dir"

    # Copy input data
    cp $pwd/02_wgdi/genepairs_info.tsv ./
    cp $pwd/02_wgdi/block.tsv ./

    echo -e "[$(date +"%T")] Calculating slope..."
    cat genepairs_info.tsv | awk -F"\t" '{OFS="\t"}NR>1{x=($5-$4)/2+$4;y=($9-$8)/2+$8;slope=(y/x); print $0"\t"x"\t"y"\t"slope}' | sed '1iBlock\tgeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00\tvalue_x\tvalue_y\tslope' > slope_genepairs_info.tsv

    echo -e "[$(date +"%T")] Running Tukey filtering on ks_NG86..."
    Rscript $pwd/bin/Filtering_Tukey.r ./ genepairs_info.tsv ${iqr} ks_NG86
    echo -e "[$(date +"%T")] Running Tukey filtering on slope..."
    Rscript $pwd/bin/Filtering_Tukey.r ./ slope_genepairs_info.tsv ${iqr} slope

    # Process ks_NG86 filtering results
    cat ks_NG86-${iqr}.genepairs_info.tsv | awk 'NR>1{if($14!="NA")print $0}' | cut -f 1-13 > aa
    grep 'Block' ks_NG86-${iqr}.genepairs_info.tsv | cat - aa | cut -f 1-13 > filtered_ks_NG86-${iqr}.genepairs_info.tsv
    rm aa

    # Process slope filtering results
    cat slope-${iqr}.genepairs_info.tsv | awk 'NR>1{if($17!="NA")print $0}' | cut -f 1-13 > aa
    grep 'Block' slope-${iqr}.genepairs_info.tsv | cat - aa | cut -f 1-13 > filtered_slope-${iqr}.genepairs_info.tsv
    rm aa slope_genepairs_info.tsv

    # Combine both filters (keep only blocks present in both)
    cat filtered_ks_NG86-${iqr}.genepairs_info.tsv filtered_slope-${iqr}.genepairs_info.tsv | grep -v 'Block' | sort | uniq -c | awk '{if($1=="2")print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sort -k1,1n | sed '1iBlock\tgeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00' > raw_filtered_ks_NG86.genepairs_info.tsv

    # Get kept block IDs from each filter
    awk 'NR>1{print $1}' filtered_ks_NG86-${iqr}.genepairs_info.tsv | sort | uniq | sort -k1,1n > keep.ksblock.id
    awk 'NR>1{print $1}' filtered_slope-${iqr}.genepairs_info.tsv | sort | uniq | sort -k1,1n > keep.positionblock.id

    # Retain blocks that passed both filters
    cat keep.ksblock.id keep.positionblock.id | sort | uniq -c | awk '{if($1==2){print $2}}' | sed '1i Block' | sort -k1,1 | awk 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a)print $0}' - raw_filtered_ks_NG86.genepairs_info.tsv > filtered.ks-slope-bioplot_genepairs_info.tsv
    rm keep.ksblock.id keep.positionblock.id

    # Final allele pairs
    cat filtered.ks-slope-bioplot_genepairs_info.tsv | awk 'NR>1{print $2"\t"$6}' > allele_pairs_lasted.txt

    # Plotting
    mkdir -p allele_plot && cd allele_plot
    for i in `cat $pwd/00_data/$chrpairs | cut -f 3`; do
        awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$1]}' $pwd/01_genetribe/${i}/${i}A.bed $work_dir/allele_pairs_lasted.txt | awk '{if(NF=="8")print $0}' > ${i}.aa
        awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$2]}' $pwd/01_genetribe/${i}/${i}B.bed ${i}.aa | awk -va=${i} '{print a"\t"$4"\t"$5"\t"$10"\t"$11}' > ${i}.coord.allele
    done
    cat *.coord.allele > pairs.coord.allele
    rm *.aa

    Rscript $pwd/bin/plot_pairs_allele.r ./ pairs.coord.allele

    # Cleanup
    cd $work_dir
    rm raw_filtered_ks_NG86.genepairs_info.tsv ks_NG86-${iqr}.genepairs_info.tsv slope-${iqr}.genepairs_info.tsv

    echo -e "[$(date +"%T")] Completed IQR multiple: $iqr (output in $work_dir)\n"
done

echo -e "\nAll IQR multiples processed. Results are in separate directories:"
for iqr in "${iqr_values[@]}"; do
    echo -e "$pwd/02_wgdi-Tukey_${iqr}"
done
