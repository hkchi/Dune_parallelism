# parallel_regions

* CLR statistic

Generate allele frequency input file and calculate genetic position of SNPs on each chromosome
```
while read chr
do
	bcftools view PetFal.snps_bi.vcf.gz -S PetFal_Dune.sample_list.txt -r $chr -O v | perl vcf2AF.perennials.pl --snp PetFal_Dune.$chr > PetFal_Dune.$chr.freq.txt
	perl SF2_map.pl PetFal_Dune.$chr.snps geneticmap.txt > PetFal_Dune.$chr.map
done < chromosome.list
```

Generate frequency spectrum across all chromosomes and compute the empirical frequency spectrum
```
cat PetFal_Dune.*.freq.txt | awk 'NR==1||$1!="position"{print}' > PetFal_Dune.CombinedFreq.txt
SweepFinder2 -f PetFal_Dune.CombinedFreq.txt PetFal_Dune.Spect.txt
```

Run SweepFinder2 (no B-value map provided)
```
while read chr
do
	SweepFinder2 -lrg 2000 PetFal_Dune.$chr.freq.txt PetFal_Dune.Spect.txt PetFal_Dune.$chr.map PetFal_Dune.$chr.SF2.txt
done < chromosome.list
```

Combine the output
```
Rscript SF2.R
```

* Reduction of diversity (ROD)

Calculate pi using VCFtools
```
window_size=20000
step=2000
vcftools --gzvcf PetFal.snps_bi.vcf.gz --keep PetFal_nonDune.sample_list.txt --window-pi $window_size --window-pi-step $step --out PetFal.nonDune
vcftools --gzvcf PetFal.snps_bi.vcf.gz --keep PetFal_Dune.sample_list.txt --window-pi $window_size --window-pi-step $step --out PetFal.Dune
```

Calculate pi_ratio
```
Rscript pi_ratio.R
```

* Genetic differentiation

Calculate Fst using VCFtools
```
window_size=20000
step=2000
vcftools --gzvcf PetFal.snps_bi.vcf.gz --weir-fst-pop PetFal_nonDune.sample_list.txt --weir-fst-pop PetFal_Dune.sample_list.txt --fst-window-size $window_size --fst-window-step $step --out PetFal
```

Format the output
```
Rscript Fst.R
```

* Regions of parallel evolution and enrichment analyses

Overlap three statistics for putative selected regions
```
# Define selected regions
perl merge_CLR_region.pl PetFal_Dune.CLR.txt PetFal_Dune
perl merge_pi_ratio_region.pl PetFal.pi_ratio.txt PetFal
perl merge_fst_region.pl PetFal.fst.txt PetFal
perl merge_CLR_region.pl PetFal2_Dune2.CLR.txt PetFal2_Dune2
perl merge_pi_ratio_region.pl PetFal2.pi_ratio.txt PetFal2
perl merge_fst_region.pl PetFal2.fst.txt PetFal2
# Find overlapping selected regions
grep -v "^#" PetFal_Dune.CLR.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal_Dune.CLR.bed
grep -v "^#" PetFal.pi_ratio.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal.pi_ratio.bed
grep -v "^#" PetFal.fst.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal.fst.bed
bedtools intersect -a PetFal.pi_ratio.bed -b PetFal.fst.bed | bedtools intersect -a PetFal_Dune.CLR.bed -b - -wa | uniq > PetFal_Dune.selected_region.bed
grep -v "^#" PetFal2_Dune2.CLR.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal2_Dune2.CLR.bed
grep -v "^#" PetFal2.pi_ratio.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal2.pi_ratio.bed
grep -v "^#" PetFal2.fst.selected_region | cut -f 1,3,4 | awk -v OFS='\t' '{if($1<10){sub(/^/,"0",$1)};sub(/^/,"Ha412HOChr",$1);print $0}' > PetFal2.fst.bed
bedtools intersect -a PetFal2.pi_ratio.bed -b PetFal2.fst.bed | bedtools intersect -a PetFal2_Dune2.CLR.bed -b - -wa | uniq > PetFal2_Dune2.selected_region.bed
# Parallel regions
bedtools intersect -a PetFal_Dune.selected_region.bed -b PetFal2_Dune2.selected_region.bed > parallel.selected_region.bed
```

LD pruning of the selected regions
```
Rscript LD_prune.R
```

Enrichment analysis of the selected regions
```
Rscript enrichment_cM.R
```

* Phylogeny of parallel regions

Phylogenetic analyses for parallel regions using Twisst
```
# Phase the VCF with Beagle 4
java -Xmx180G -jar beagle.18May20.d20.jar \
	gt=allPET.bi.vcf \
	out=allPET.bi.beagle \
	impute=true \
	window=1 \
	overlap=0.1 \
	nthreads=20
bcftools index allPET.bi.beagle.vcf.gz
# ML trees for genomic windows
mkdir ./iqtree/
region_file="phylo_region.bed"
win=20000
while read chr_i start_i end_i
do
	region_name="${chr_i}_${start_i}_${end_i}"
	cd ./iqtree/
	if [ ! -s ${chr_i}.vcf.gz ]
	then
		bcftools view -S ../sample_list -r ${chr_i} ../allPET.bi.beagle.vcf.gz -O z -o ${chr_i}.vcf.gz
		bcftools index ${chr_i}.vcf.gz
	fi
	# Genomic window in the region
	grep ${chr_i} ../HA412.chr_len.txt | bedtools makewindows -g - -w $win > tmp.${chr_i}.win${win}.bed
	echo -e "${chr_i}\t${start_i}\t${end_i}" | bedtools intersect -a tmp.${chr_i}.win${win}.bed -b - -f 0.5 -wa > tmp.$region_name.win${win}.bed
	# Run iqtree for each window
	while read chr start end
	do
		bcftools view ${chr_i}.vcf.gz -r ${chr}:${start}-${end} -O v | perl ../vcf2fasta.twisst.pl tmp.${chr}_${start}_${end}
	done < tmp.$region_name.win${win}.bed
	if [ -s $region_name.win$win.win_pos.txt ]
	then
		rm $region_name.win$win.win_pos.txt
	fi
	k=0
	while read chr start end
	do
		if [ -s tmp.${chr}_${start}_${end}.snps.fasta ]
		then
			((k=k+1))
			echo -e "$k\t${chr}\t${start}\t${end}" >> $region_name.win$win.win_pos.txt
			iqtree -s tmp.${chr}_${start}_${end}.snps.fasta -o DIV_1956,GRO_2043 -st DNA -m MFP+ASC -nt 20 -redo -pre tmp.$region_name.$k
		fi
	done < tmp.$region_name.win${win}.bed
	# Gather the trees
	mv $region_name.win$win.win_pos.txt ../
	myfile=""
	for k in $(seq 1 $(wc -l ../$region_name.win$win.win_pos.txt | cut -d " " -f1))
	do
		myfile="$myfile tmp.$region_name.$k.treefile"
	done
	cat $myfile | gzip > ../$region_name.win$win.trees.gz
	#rm tmp.*
	cd ../
	# Twisst - Run 1 (modes of sources: habitat, geography, others)
	python twisst/twisst.py \
		-t $region_name.win$win.trees.gz \
		--inputTopos mode1.topos.txt \
		-w mode1.$region_name.win$win.weights.txt \
		--groupsFile dune_phylogeny.list \
		-g GSD_dune -g GSD_nondune -g MON_dune -g MON_nondune \
		--method complete
	# Twisst - Run 2 (sources within species: within species diversity, introgression/ancient variation, subspecies introgression)
	python twisst/twisst.py \
		-t $region_name.win$win.trees.gz \
		--inputTopos mode2_1.topos.txt \
		-w mode2_1.$region_name.win$win.weights.txt \
		--groupsFile dune_phylogeny.list \
		-g GSD_dune -g GSD_nondune -g petpet -g outgroup \
		--method complete
	python twisst/twisst.py \
		-t $region_name.win$win.trees.gz \
		--inputTopos mode2_2.topos.txt \
		-w mode2_2.$region_name.win$win.weights.txt \
		--groupsFile dune_phylogeny.list \
		-g MON_dune -g MON_nondune -g petpet -g outgroup \
		--method complete
done < $region_file
```
