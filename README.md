#work directory
/well/ckb/users/ahd580/PGS
#software 
/well/ckb/users/ahd580/tools/gctb_2.02_Linux/gctb #SbayeasR
/well/ckb/users/ahd580/tools/PRScs/PRScs.py #PRS-CS
/well/ckb/users/ahd580/tools/PRS2.3.3/ #PRsice2.3.3
qrsh -q short.qc -P ckb.prjc -pe shmem 2

####dataset
## CKB bfile /well/ckb/shared/ckb_bfile
qsub /well/ckb/shared/ckb_bfile/ckb_bfile.sh #convert begn to bfile
qsub /well/ckb/shared/ckb_bfile/ckb_filter.sh #filter out duplicated SNPs && SNPs with low quality i.e. maf <0.005, info score <0.3

awk '{if ($5<$6) {print $1"_"$4"_"$5"_"$6,$2} else print $1"_"$4"_"$6"_"$5,$2}' \
/well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chr*_v3.bim \
>/well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chrall.snps #snpid chr_bp_a1_a2 rsid

awk '{if ($5<$6) {print $1"_"$4"_"$5"_"$6,$2} else print $1"_"$4"_"$6"_"$5,$2}' \
/well/ckb/shared/ckb_bfile/newCHR*.bim \
>/well/ckb/shared/ckb_bfile/CHRall.snps #snpid chr_bp_a1_a2 rsid

join -1 1 -2 1 -o 1.2,2.2 <(sort -k1,1 -u /well/ckb/shared/ckb_bfile/CHRall.snps) \
<(sort -k2,2 -u /well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chrall.snps|sort -k1,1 -u ) |sort -k1,1 -u \
>/well/ckb/shared/ckb_bfile/CHRall.snps.commn.ukbsnpid #obtain snps in common for snp name update

qsub /well/ckb/shared/ckb_bfile/ckb_bfile.newid.sh ###update ckb snpid to ukb marker id (i.e. rsid)
##ukb reference /well/ckb/users/ahd580/tools/ukb15kref

for chrom in $(seq 1 22 );do
awk '{print $2}'  /well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chr${chrom}_v3.bim|sort|uniq -d \
>/well/ckb/users/ahd580/tools/ukb15kref/${chrom}.snp.exl

/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chr${chrom}_v3 \
--maf 0.01 \
--geno 0.05 \
--exclude /well/ckb/users/ahd580/tools/ukb15kref/${chrom}.snp.exl \
--make-bed \
--out /well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr${chrom}
 done #filter: duplicated snps, maf0.01 call rate 0.95

##phenotype 
/well/ckb/users/ahd580/PGS/zammy/ckb_men_pheno_gsid.txt #n=41255
/well/ckb/users/ahd580/PGS/zammy/ckb_women_pheno_gsid.txt #n=54607
 # WHR_RINT were obtained by pooling sex and RC stratfied linear regression residuals of WHR adjusted for age and age2.
 #training data
 cat <(head -n 1 /well/ckb/users/ahd580/PGS/zammy/ckb_men_pheno_gsid.txt ) \
 <( grep -Fwf <(awk '{print $1}' /well/ckb-share/ckb_v1.0/unrelated/national.king.cutoff.in.id) \
 /well/ckb/users/ahd580/PGS/zammy/ckb_men_pheno_gsid.txt) \
 >/well/ckb/users/ahd580/PGS/zammy/unrel.ckb_men_pheno_gsid.txt ## unrelated sub population were used for training n=29499

cat <(head -n 1 /well/ckb/users/ahd580/PGS/zammy/ckb_women_pheno_gsid.txt ) \
 <( grep -Fwf <(awk '{print $1}' /well/ckb-share/ckb_v1.0/unrelated/national.king.cutoff.in.id) \
 /well/ckb/users/ahd580/PGS/zammy/ckb_women_pheno_gsid.txt) \
 >/well/ckb/users/ahd580/PGS/zammy/unrel.ckb_women_pheno_gsid.txt  ## unrelated sub population were used for training n=42295

#test data
 grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/SHR_c.txt) \
 /well/ckb/users/ahd580/PGS/zammy/ckb_men_pheno_gsid.txt \
 >/well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_men_pheno_gsid.txt ## unrelated sub population were used for training n=7316

grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/SHR_c.txt) \
 /well/ckb/users/ahd580/PGS/zammy/ckb_women_pheno_gsid.txt \
 >/well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_women_pheno_gsid.txt  ## unrelated sub population were used for training n=8212


####-----------------------####
#            PRSice2.3.3      #
####-----------------------####

match with ckb snp markers and add column of ckb snp marker name


awk '{if ($5<$6) {print $1"_"$4"_"$5"_"$6,$2} else print $1"_"$4"_"$6"_"$5,$2}' \
/well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chr*_v3.bim \
>/well/ckb/users/ahd580/tools/ukb15kref/random15k_ref.chrall.snps #snpid chr_bp_a1_a2 rsid

awk '{if ($6<$7) {print $4"_"$5"_"$6"_"$7,$0} else print $4"_"$5"_"$7"_"$6,$0}' \
/well/ckb/users/ahd580/PGS/zammy/WHR_Women_SNP.bgen \
>/well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen

awk '{if ($6<$7) {print $4"_"$5"_"$6"_"$7,$0} else print $4"_"$5"_"$7"_"$6,$0}' \
/well/ckb/users/ahd580/PGS/zammy/WHR_Men_SNP.bgen \
>/well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen

grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen) \
<(awk '$8>0.3&&$7>0.005&&$7<0.995 {print $1,$2}' /well/ckb-share/ckb_v1.0/info.table|sort -k1,1 -u ) \
|sort -k1,1 -u \
>/well/ckb/users/ahd580/PGS/zammy/snps.commn.ckbsnpid #obtain snps in common for snp name update

cat <(echo 'rsid' $(head -n 1 /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen)) \
<(join -1 1 -2 1 -o 2.2,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13 \
<(sort -k1,1 -u /well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen) \
<(awk '$8>0.3&&$7>0.005&&$7<0.995 {print $1,$2}' /well/ckb-share/ckb_v1.0/info.table|sort -k1,1 -u ) \
|sort -k1,1 -u) \
> /well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen.txt

cat <(echo 'rsid' $(head -n 1 /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen)) \
<(join -1 1 -2 1 -o 2.2,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13 \
<(sort -k1,1 -u /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen) \
<(awk '$8>0.3&&$7>0.005&&$7<0.995 {print $1,$2}' /well/ckb-share/ckb_v1.0/info.table|sort -k1,1 -u ) \
|sort -k1,1 -u) \
> /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.txt


# creat soft link for the ld bfiles
for file in 'bed' 'bim' 'fam';do
for i in $(seq 1 22);do
chrom=$(printf %02g $i)
ln -s /well/ckb/users/bcd680/data/imputed_ld_bed_v1.1/$chrom.$file \
/well/ckb/shared/ckbld_sl/$i.$file;
done 
done
#--clump-kb 500 --clump-r2 0.05 --ld /well/ckb/users/bcd680/data/imputed_ld_bed_v1.1/ --proxy 0.8

wellPRSice=/well/ckb/users/ahd580/tools/PRS2.3.3
wd=/well/ckb/users/ahd580/PGS/zammy
phenofile=$1
gender=$2
base_gwas=$wd/WHR_${gender}_SNP.bgen.txt
plower=$(sort -k14,14g $base_gwas|awk 'NR==2{print $14}')
puper=$3
pinterv=$4
targetbgen=/well/ckb-share/ckb_v1.0/impute_qced
/well/ckb-share/local/bin/Rscript $wellPRSice/PRSice.R --prsice $wellPRSice/PRSice_linux --dir $wd \
--binary-target F \
--base $base_gwas \
--snp rsid \
--A1 a1 \
--A2 a2 \
--beta \
--pvalue P_value \
--chr CHR \
--bp BP \
--stat beta \
--pheno-file $wd/$phenofile \
--cov $wd/$phenofile \
--pheno-col WHR_RINT \
--cov-col @PC[1-10],array \
--cov-factor array \
--type bgen \
--target $targetbgen/#,$targetbgen/23.sample \
--thread 2 \
--lower $plower  \
--interval $pinterv \
--upper $puper \
--clump-kb 500 \
--clump-r2 0.05 \
--ld /well/ckb/shared/ckbld_sl/# \
--proxy 0.8 \
--base-info INFO:0.3 \
--maf 0.005 \
--info 0.3 \
--out $wd/PRsice/PRSice.$phenofile

#--bgen file format

for SEX in women men ;do
PHENOFILE=unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.sh $PHENOFILE $SEX 1 1e-5
done

for SEX in women men ;do
PHENOFILE=ex.unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.testing.sh $PHENOFILE $SEX $(awk 'NR==2 {print $3}' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.unrel.ckb_${SEX}_pheno_gsid.txt.summary) 1
done
#--bed file format

#!/bin/bash
#files name ckbld.sh
#$ -cwd -V -N CKBLD
#$ -q short.qc
#$ -P ckb.prjc
#$ -t 1-22
chrom=$(printf %01g $SGE_TASK_ID)
/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/shared/ckb_bfile/chr${chrom} \
--keep /well/ckb/shared/ckbld_sl/22.fam \
--make-bed --out /well/ckb/shared/ckb_bfile/ld/${chrom}

for SEX in women men ;do
PHENOFILE=unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.bfile.sh $PHENOFILE $SEX 1 1e-5
done

for SEX in women men ;do
PHENOFILE=ex.unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.testing.bfile.sh $PHENOFILE $SEX $(awk 'NR==2 {print $3}' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.unrel.ckb_${SEX}_pheno_gsid.txt.bed.summary) 1
done

##GRS
for SEX in women men ;do
PHENOFILE=unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.GRS.sh $PHENOFILE $SEX 5e-8 1e-100
done

for SEX in women men ;do
PHENOFILE=ex.unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.GRS.testing.sh $PHENOFILE $SEX $(awk 'NR==2 {print $3}' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.GRS.unrel.ckb_${SEX}_pheno_gsid.txt.summary) 1
done
#--bed file format
for SEX in women men ;do
PHENOFILE=unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.GRS.bfile.sh $PHENOFILE $SEX 5e-8 1e-100
done

for SEX in women men ;do
PHENOFILE=ex.unrel.ckb_${SEX}_pheno_gsid.txt
qsub PRsice2.3.3.GRS.testing.bfile.sh $PHENOFILE $SEX \
$(awk '$3<5e-8' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.GRS.unrel.ckb_${SEX}_pheno_gsid.txt.bed.prsice|sort -rk4,4g|tail -n 1|awk '{print $3}') 1
done

echo 'Data Sex Threshold PRS.R2 Coefficient Standard.Error P Num_SNP' > /well/ckb/users/ahd580/PGS/zammy/PRsice/PRsice2.3.3.results.txt
for grs in '' ;do
for sub in 'unrel.' ;do
for SEX in men women;do
echo $(awk 'NR==2 {print  "'$grs'","'$sub'","'$SEX'",$3,$8,$9,$10,$4,$11}' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.${grs}${sub}ckb_${SEX}_pheno_gsid.txt.bed.summary) >>/well/ckb/users/ahd580/PGS/zammy/PRsice/PRsice2.3.3.results.txt
done
done
done

for grs in 'GRS.';do
for sub in 'unrel.' ;do
for SEX in men women;do
echo $(awk '$3<5e-8' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.${grs}${sub}ckb_${SEX}_pheno_gsid.txt.bed.prsice|sort -rk4,4g|tail -n 1| \
awk '{print  "'$grs'","'$sub'","'$SEX'",$3,$6,$7,$5,$4,$8}' ) >>/well/ckb/users/ahd580/PGS/zammy/PRsice/PRsice2.3.3.results.txt
done
done
done


for grs in '' 'GRS.';do
for sub in  'ex.unrel.';do
for SEX in men  women;do
echo $(awk 'NR==2 {print  "'$grs'","'$sub'","'$SEX'",$3,$6,$7,$5,$4,$8}' /well/ckb/users/ahd580/PGS/zammy/PRsice/PRSice.${grs}${sub}ckb_${SEX}_pheno_gsid.txt.bed.prsice) >>/well/ckb/users/ahd580/PGS/zammy/PRsice/PRsice2.3.3.results.txt
done
done
done

###GRS using list of SNPs and beta from zammy

join -1 1 -2 1 -o 2.2,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.7,2.8 <(awk '{if ($4<$5) {print $2"_"$3"_"$4"_"$5,$0} else print $2"_"$3"_"$5"_"$4,$0}' /well/ckb/users/tpp788/scores/snp_lists/UKB_WHR_extract.txt|sed 's/Chromosome_Position_A1_A2/loc_aa/'|sort -k1,1) <(sort -k1,1 /well/ckb-share/ckb_v1.0/info.table)|sort -k4g,5g|sed 's/Chromosome/CHR/'|sed 's/Position/BP/'|sed 's/info/INFO/'|sed 's/BETA_Men/BETA_men/'|sed 's/BETA_Women/BETA_women/'|sed 's/NA/1/'| awk '{if ($1=="rsid") {print $0,"P_men P_women"} else print $0,1,1}' >/well/ckb/users/ahd580/PGS/zammy/UKB_WHR_extract.rsid.txt

for SEX in women men ;do
for PHENOFILE in ex.unrel.ckb_${SEX}_pheno_gsid.txt unrel.ckb_${SEX}_pheno_gsid.txt;do
qsub PRsice2.3.3.GRS386.sh $PHENOFILE $SEX 1 1
done
done

 ####-------------------------####
 #               SBayesR         #
 ####-------------------------####


 #reformat the gwas summary stats as .ma format
#SNP A1 A2 freq b se p N 
cat <(echo 'SNP A1 A2 freq b se p N') \
<(awk 'NR>1{print $2,$6,$7,$8,$10,$11,$12,203124}' /well/ckb/users/ahd580/PGS/zammy/WHR_Men_SNP.bgen) \
>/well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.ma

cat <(echo 'SNP A1 A2 freq b se p N') \
<(awk '{print $2,$6,$7,$8,$10,$11,$12,240145}' /well/ckb/users/ahd580/PGS/zammy/WHR_Women_SNP.bgen) \
>/well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen.ma

##
chrom=$(printf %01g $SGE_TASK_ID)
wellSBayesR=/well/ckb/users/ahd580/tools/gctb_2.02_Linux
wd=/well/ckb/users/ahd580/PGS/zammy
gender=$1
out=$2
gwassum=$wd/WHR_${gender}_SNP.bgen.ma
targetbfile=/well/ckb/shared/ckb_bfile

$wellSBayesR/gctb --sbayes R \
--ldm $wellSBayesR/ukb_50k_bigset_2.8M/ukb50k_shrunk_chr${chrom}_mafpt01.ldm.sparse \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--gwas-summary $gwassum \
--chain-length 10000 \
--burn-in 2000 \
--out-freq 10 \
--out $wd/$out


for SEX in women men ;do
qsub SBayesR.sh $SEX
done


#!/bin/bash
#files name SBayesR.plinkscore.sh
#$ -cwd -V -N TrainningSBayesR.score
#$ -q long.qc
#$ -P ckb.prjc
#$ -t 1-22
chrom=$(printf %01g $SGE_TASK_ID)
gender=$1

/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/shared/ckb_bfile/chr${chrom} \
--score /well/ckb/users/ahd580/PGS/zammy/SBayesR/${gender}.chr${chrom}.snpRes 2 5 8 sum \
--out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_${gender}.chr${chrom}.score

for SEX in women men ;do
qsub SBayesR.plinkscore.sh $SEX 
done

##R score
for (sub in c('unrel', 'ex.unrel')){
    for (sex in c('women','men')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_',sex,'.chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in rep(c(1000,2000,4000,6000,7000,nrow(pheno),10000,15000,20000,25000),each=5)){
    fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
print (paste0(sub,' ',n,' ',sex,' ',summary(fit)$r.squared))
print (summary(fit)$coefficients[2,c(1,2,4)])
}
}
}

###----- using ckb ld references -------------###

##Genetc map files
# download ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CHB_omni_recombination_20130507.tar
# /well/ckb/users/ahd580/tools/gctb_2.02_Linux/
for chrom in $(seq 1 22);do
 join -1 2 -2 1 -o 1.1,2.1,2.2 <(awk '{print $2,$4}' /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.bim|sort -k2,2) \
 <(awk '$4==0 {print $1,$3}' /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/CHB-${chrom}-final.txt|sort -k1,1) \
>/well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/CHB-${chrom}-final.genetic.map
done

for chrom in $(seq 1 22);do
v=$(grep -Fwf <(awk '{print $2 }' /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.bim) /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/CHB-${chrom}-final.genetic.map|wc -l|awk '{print $1}')
c=$(expr $v / 5000 - 1)
for i in $(seq 1 $c); do 
echo "$(( i * 5000 -4999))-$(( i * 5000 ))" \
>>/well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange
done
echo "$(( c * 5000 + 1 ))-$v" \
>>/well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange
done

for chrom in $(seq 1 22); do
job=$(wc -l /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
ckbrange=$(awk -v n=$j 'NR==n' /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange)
printf '##SBayesR_CKB_LDM\n#!/bin/bash\n''#SBayesR_ckb_ldm_chr'$chrom'_'$j'.sh \n''#$ -cwd -V -N SBayesR_ckb_ldm_chr'$chrom'_'$j'\n''#$ -q short.qc\n''#$ -P ckb.prjc\n''#$ -pe shmem 1\n''/well/ckb-share/local/bin/gctb --bfile /well/ckb/shared/ckb_bfile/ld/CKB_hm3'$chrom' --make-shrunk-ldm --snp '$ckbrange' --gen-map /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/CHB-'$chrom'-final.genetic.map --ne 13038.18 --genmap-n 108 --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr'$chrom'\n' \
>/well/ckb/users/ahd580/PGS/zammy/script/SBayesR_ckb_ldm_chr${chrom}_${j}.sh
done
done

#------example---------#
##SBayesR_CKB_LDM
#!/bin/bash
#SBayesR_ckb_ldm_chr22_1.sh
#$ -cwd -V -N SBayesR_ckb_ldm_chr22_1
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 1
/well/ckb-share/local/bin/gctb --bfile /well/ckb/shared/ckb_bfile/ld/CKB_hm322 --make-shrunk-ldm --snp 1-9061 --gen-map /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/CHB-22-final.genetic.map --ne 13038.18 --genmap-n 108 --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr22

for chrom in $(seq 1 22); do
job=$(wc -l /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
qsub /well/ckb/users/ahd580/PGS/zammy/script/SBayesR_ckb_ldm_chr${chrom}_${j}.sh
done
done

for chrom in $(seq 1 22); do
rm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${chrom}.mldmlist
job=$(wc -l /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
ckbrange=$(awk -v n=$j 'NR==n' /well/ckb/shared/ckb_bfile/ld/CKB_hm3${chrom}.ldm.snprange)
echo '/well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr'${chrom}'.snp'${ckbrange}'.ldm.shrunk' \
>>/well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${chrom}.mldmlist
done
done

for j in $(seq 1 22); do
/well/ckb/users/ahd580/tools/gctb_2.02_Linux/gctb --mldm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${j}.mldmlist \
--make-full-ldm --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${j}
done

for chrom in $(seq 1 22); do
qsub SBayesR_ckb_ldm_join.sh $chrom
done



for j in $(seq 1 22); do
/well/ckb-share/local/bin/gctb \
--ldm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${j}.ldm.full --make-sparse-ldm --chisq 0 \
--out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ckb_ldm_chr${j}.ldm.full.sparse
done

qsub SBayesR_ckb_ldm_join_sparse.sh

for SEX in women men ;do
qsub SBayesR_ckbld.sh $SEX
done

for SEX in women men ;do
qsub SBayesR.plinkscore.ckbld.sh $SEX 
done

for (sub in c('unrel', 'ex.unrel')){
    for (sex in c('women','men')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/SBayesR/ckbld_ckb_',sex,'.chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in nrow(pheno)){
    fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
print (paste0(sub,' ',n,' ',sex,' ',summary(fit)$r.squared))
print (summary(fit)$coefficients[2,c(1,2,4)])
}
}
}

###----- using ukb 15k ld references -------------###

##Genetc map files
# download ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar
# /well/ckb/users/ahd580/tools/gctb_2.02_Linux/
for chrom in $(seq 1 22);do
 join -1 2 -2 1 -o 1.1,2.1,2.2 <(awk '{print $2,$4}' /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.bim|sort -k2,2) \
 <(awk '$4==0 {print $1,$3}' /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/CEU-${chrom}-final.txt|sort -k1,1) \
>/well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/CEU-${chrom}-final.genetic.map
done

 rm /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3*.ldm.snprange

for chrom in $(seq 1 22);do
v=$(grep -Fwf <(awk '{print $2 }' /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.bim) /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/CEU-${chrom}-final.genetic.map|wc -l|awk '{print $1}')
c=$(expr $v / 5000 - 1)
for i in $(seq 1 $c); do 
echo "$(( i * 5000 -4999))-$(( i * 5000 ))" \
>>/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange
done
echo "$(( c * 5000 + 1 ))-$v" \
>>/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange
done


for chrom in $(seq 1 22); do
job=$(wc -l /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
ukbrange=$(awk -v n=$j 'NR==n' /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange)
printf '##SBayesR_UKB_LDM\n#!/bin/bash\n''#SBayesR_ukb_ldm_chr'$chrom'_'$j'.sh \n''#$ -cwd -V -N SBayesR_ukb_ldm_chr'$chrom'_'$j'\n''#$ -q short.qc\n''#$ -P ckb.prjc\n''#$ -pe shmem 1\n''/well/ckb-share/local/bin/gctb --bfile /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3'$chrom' --make-shrunk-ldm --snp '$ukbrange' --gen-map /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/CEU-'$chrom'-final.genetic.map --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr'$chrom'\n' \
>/well/ckb/users/ahd580/PGS/zammy/script/SBayesR_ukb_ldm_chr${chrom}_${j}.sh
done
done

for chrom in $(seq 1 22); do
job=$(wc -l /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
qsub /well/ckb/users/ahd580/PGS/zammy/script/SBayesR_ukb_ldm_chr${chrom}_${j}.sh
done
done

#------example---------#
##SBayesR_UKB_LDM
#!/bin/bash
#SBayesR_ukb_ldm_chr22.sh
#$ -cwd -V -N SBayesR_ukb_ldm_chr22
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 2
#$ -t 1-1
j=$(printf %01g $SGE_TASK_ID)
ukbrange22=$(awk -v n=$j 'NR==n' /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm322.ldm.snprange)
/well/ckb-share/local/bin/gctb --bfile /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm322 --make-shrunk-ldm --snp $ukbrange22 --gen-map /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/CEU-22-final.genetic.map --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr22

for chrom in $(seq 1 22); do
qsub /well/ckb/users/ahd580/PGS/zammy/script/SBayesR_ukb_ldm_chr${chrom}.sh
done


for chrom in $(seq 1 22); do
rm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${chrom}.mldmlist
job=$(wc -l /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange|cut -d' ' -f1)
for j in $(seq 1 $job);do
ukbrange=$(awk -v n=$j 'NR==n' /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3${chrom}.ldm.snprange)
echo '/well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr'${chrom}'.snp'${ukbrange}'.ldm.shrunk' \
>>/well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${chrom}.mldmlist
done
done

#!/bin/bash
#SBayesR_ukb_ldm_join.sh
#$ -cwd -V -N SBayesR_ukb_ldm_join
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 2
#$ -t 1-22
j=$(printf %01g $SGE_TASK_ID)
/well/ckb-share/local/bin/gctb --mldm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${j}.mldmlist \
--make-full-ldm --out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${j}

qsub SBayesR_ukb_ldm_join.sh

#!/bin/bash
#SBayesR_ukb_ldm_join_sparse.sh
#$ -cwd -V -N SBayesR_ukb_ldm_join_sparse
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 2
#$ -t 1-22
j=$(printf %01g $SGE_TASK_ID)
/well/ckb-share/local/bin/gctb \
--ldm /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${j}.ldm.full --make-sparse-ldm --chisq 0 \
--out /well/ckb/users/ahd580/PGS/zammy/SBayesR/ukb_ldm_chr${j}.ldm.full.sparse

qsub SBayesR_ukb_ldm_join_sparse.sh

##PRScs
#summary file format # SNP          A1   A2   BETA      P
cat <(echo 'SNP A1 A2 BETA P') \
<(awk 'NR>1{print $2,$6,$7,$10,$12}' /well/ckb/users/ahd580/PGS/zammy/WHR_Men_SNP.bgen) \
>/well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.prscs.txt
## creat training/testing bfile seperately
grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/unrel.ckb_men_pheno_gsid.txt) /well/ckb/shared/ckb_bfile/chr1.fam |awk '{print $1,$2}' >/well/ckb/users/ahd580/PGS/zammy/PRScs/unrel.ckb_men.sample.txt
grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/unrel.ckb_women_pheno_gsid.txt) /well/ckb/shared/ckb_bfile/chr1.fam |awk '{print $1,$2}' >/well/ckb/users/ahd580/PGS/zammy/PRScs/unrel.ckb_women.sample.txt
 grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_men_pheno_gsid.txt) /well/ckb/shared/ckb_bfile/chr1.fam |awk '{print $1,$2}' >/well/ckb/users/ahd580/PGS/zammy/PRScs/ex.unrel.ckb_men.sample.txt
grep -Fwf <(awk '{print $1}' /well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_women_pheno_gsid.txt) /well/ckb/shared/ckb_bfile/chr1.fam |awk '{print $1,$2}' >/well/ckb/users/ahd580/PGS/zammy/PRScs/ex.unrel.ckb_women.sample.txt

for subdata in 'ex.unrel' 'unrel';do
for sex in women men;do
for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/shared/ckb_bfile/chr${chrom} \
--keep /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.sample.txt \
--make-bed \
--out /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.chr${chrom}
done
done
done

cat <(echo 'SNP A1 A2 BETA P') \
<(awk '{print $2,$6,$7,$10,$12}' /well/ckb/users/ahd580/PGS/zammy/WHR_Women_SNP.bgen) \
>/well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen.prscs.txt
cat <(echo 'SNP A1 A2 BETA P') \
<(awk '{print $2,$6,$7,$10,$12}' /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen) \
>/well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.prscs.txt

#!/bin/bash
#files name PRScs.sh
#$ -cwd -V -N TrainningPRScs
#$ -q long.qc
#$ -P ckb.prjc
#$ -pe shmem 6
#$ -t 1-22
chrom=$(printf %01g $SGE_TASK_ID)
sex=$1
n=$2
python /well/ckb/users/ahd580/tools/PRScs/PRScs.py \
--ref_dir=/well/ckb/users/ahd580/tools/PRScs/ldblk_1kg_eur/ \
--bim_prefix=/well/ckb/users/ahd580/PGS/zammy/PRScs/unrel.ckb_${sex}.chr${chrom} \
--sst_file=/well/ckb/users/ahd580/PGS/zammy/WHR_${sex}_SNP.bgen.prscs.txt \
--n_gwas=$n \
--out_dir=/well/ckb/users/ahd580/PGS/zammy/PRScs/$sex/


qsub PRScs.sh men 203124

qsub PRScs.sh women 240145

qsub PRScs.EAS.sh men 203124

qsub PRScs.EAS.sh women 240145


#!/bin/bash
#files name PRScs.plinkscore.sh
#$ -cwd -V -N TrainningPRScs.score
#$ -q long.qc
#$ -P ckb.prjc
#$ -t 1-22
chrom=$(printf %01g $SGE_TASK_ID)
gender=$1

/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/shared/ckb_bfile/chr${chrom} \
--score /well/ckb/users/ahd580/PGS/zammy/PRScs/${gender}/_pst_eff_a1_b0.5_phiauto_chr${chrom}.txt 2 4 6 sum \
--out /well/ckb/users/ahd580/PGS/zammy/PRScs/${gender}/chr${chrom}.score

for SEX in women men ;do
qsub PRScs.plinkscore.sh $SEX 
done

for SEX in women men ;do
qsub PRScs.plinkscore.EAS.sh $SEX 
done


##R score
result=as.data.frame(array(dim = c(1,7),NA))
for (sub in c('unrel', 'ex.unrel')){
    for (sex in c('women','men')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in rep(c(1000,2000,4000,6000,7000,nrow(pheno),10000,15000,20000,25000),each=5)){
fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
result=rbind(result,c(sub,n,sex,summary(fit)$r.squared,summary(fit)$coefficients[2,c(1,2,4)]))
}
}
}
 result=result[-1,]    
names(result)=c('Data','N','Sex','r2','BETA','SE','P')
result=result[with(result,order(Data,Sex,N)),]

##investigate variation of the PRS explained R2 explained by sample size 
result=NULL
for (sub in c('unrel')){
    for (sex in c('women')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
sample=read.table('/well/ckb-share/ckb_v1.0/sample.info',header=T)
pheno=merge(pheno,sample[,c('ID_1','RC')],by.x='FID',by.y='ID_1',all.x=T)

for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,25:46])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno)))

for (r in 1:5){
for (n in 1:5){
pheno$group<- sample(c(rep(seq(1:n),each=floor(nrow(pheno)/n)),sample(1:n,nrow(pheno)%%n,)),nrow(pheno))
for (j in 1:n){
fit=lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno[ pheno$group==j,])
fit2=lm(WHR_RINT~prs+array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno[ pheno$group==j,])
result=rbind(result,c(r,n,j,nrow(pheno[ pheno$group==j,]),sex,summary(fit2)$r.squared-summary(fit)$r.squared))
}
}
}
}
}
result=as.data.frame(result)
names(result)=c('Repeat','Group','Sample','N','Sex','r2')
result=result[with(result,order(Group,Repeat)),]


result1<-NULL
for (sub in c('unrel')){
    for (sex in c('women')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
sample=read.table('/well/ckb-share/ckb_v1.0/sample.info',header=T)
pheno=merge(pheno,sample[,c('ID_1','RC')],by.x='FID',by.y='ID_1',all.x=T)

for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,25:46])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno)))

for (r in 1:5){
for (k in 1:10){
n=floor(nrow(pheno)/845)
pheno$group<- sample(c(rep(seq(1:n),each=floor(nrow(pheno)/n)),rep(1+n,nrow(pheno)%%n)))
five_group_id<-sample(1:n,5*k)
pheno2=subset(pheno,subset=group %in% five_group_id)
fit=lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno2)
fit2=lm(WHR_RINT~prs+array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=pheno2)
result1=rbind(result1,c(r,k,'sum_5',nrow(pheno2),sex,summary(fit2)$r.squared-summary(fit)$r.squared))

for (j in 1:5){
    five_group_id2=as.data.frame(split(five_group_id, 1:5)[j])[,1]
fit=lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=subset(pheno2,subset=group %in%  five_group_id2))
fit2=lm(WHR_RINT~prs+array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+RC,data=subset(pheno2,subset=group %in%  five_group_id2))
result1=rbind(result1,c(r,k,j,nrow(subset(pheno2,subset=group %in%  five_group_id2)),sex,summary(fit2)$r.squared-summary(fit)$r.squared))
}
}
}
}
}
result1=as.data.frame(result1)
names(result1)=c('Repeat','Group','Subgroup','N','Sex','r2')
result1=result1[with(result1,order(Group,Repeat)),]

write.csv(result1,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result1_by_groups.csv')

result1_sum<- result1%>% mutate_at(.vars=c( "N", "r2"),.funs=as.numeric)%>% group_by(Group,Repeat,N) %>% summarise_at(.vars=c( "r2"), .funs=c(mean,sd))
result1_sum=as.data.frame(result1_sum)
names(result1_sum)[4:5]=c( 'r2_mean', 'r2_sd')

write.csv(result1,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result1_by_groups.csv')
write.csv(result1_sum,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result1_by_groups_sum.csv')
###

result2=NULL
for (sub in c('unrel')){
    for (sex in c('women')) {
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in rep(c(round(nrow(pheno)/1),round(nrow(pheno)/2),round(nrow(pheno)/3),round(nrow(pheno)/4),round(nrow(pheno)/5)),each=5)){
fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
result2=rbind(result2,c(n,sex,summary(fit)$r.squared,summary(fit)$coefficients[2,c(1,2,4)]))
}
}
}
result2=as.data.frame(result2)
names(result2)=c('N','Sex','r2','BETA','SE','P')

#summarize the mean and sd of the reuslts
library(dplyr)
result_sum<- result%>% mutate_at(.vars=c( "N", "r2"),.funs=as.numeric)%>% group_by(Repeat,Group) %>% summarise_at(.vars=c( "N", "r2"), .funs=c(mean,sd))
result_sum=as.data.frame(result_sum)
names(result_sum)[3:6]=c( 'N_mean', 'r2_mean','N_sd', 'r2_sd')

result_sum=result_sum[with(result_sum,order(Repeat,Group)),]

write.csv(result,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result_by_groups.csv')
write.csv(result_sum,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result_by_groups_sum.csv')

result2_sum<- result2%>% mutate_at(.vars=c( "N", "r2"),.funs=as.numeric)%>% group_by(N) %>% summarise_at(.vars=c( "r2"), .funs=c(mean,sd))
result2_sum=as.data.frame(result2_sum)
names(result2_sum)[2:3]=c( 'r2_mean', 'r2_sd')

write.csv(result2,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result2_by_groups.csv')
write.csv(result2_sum,'/well/ckb/users/ahd580/PGS/zammy/PRScs/result2_by_groups_sum.csv')

#barplot 
# (2) Bar plots + upper error bars.
pdf(file =paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/PRS-CS_R2_by_sample_size_sampling.pdf'),width = 24,height = 12 )
result<-read.csv('/well/ckb/users/ahd580/PGS/zammy/PRScs/result1_by_groups.csv')
result=result[with(result,order(Group,Repeat)),]
result$order=1:nrow(result)
ggplot(result, aes(order, r2)) +
  geom_bar(aes(fill =Subgroup  ), stat = "identity",
           position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = r2, ymax =r2, group =Subgroup ),
    width = 0.2, position = position_dodge(0.8)
  )+
  ylab('r2') +ylab('')+
			facet_wrap(vars(Group,Repeat), ncol=10,nrow=5,strip.position = "bottom", scales = "free_x",shrink=TRUE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
		 geom_text(aes(label = round(r2,3)), 
            position = position_dodge(width = 0.9), size =3,vjust = -0.25)+
					 geom_text(aes(label = N,y=-0.005), 
            position = position_dodge(width = 0.9),size =2.5, vjust = -0.25)
			
            dev.off()

pdf(file =paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/PRS-CS_R2_by_sample_size_split.pdf'),width = 24,height = 12 )
result<-read.csv('/well/ckb/users/ahd580/PGS/zammy/PRScs/result_by_groups.csv')
result=result[with(result,order(Group,Repeat)),]
result=result[-c(2:5),]
result$order=1:nrow(result)
ggplot(result, aes(order, r2)) +
  geom_bar(aes(fill =Sample  ), stat = "identity",
           position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = r2, ymax =r2, group =Sample ),
    width = 0.2, position = position_dodge(0.8)
  )+
  ylab('r2') +ylab('')+
			facet_wrap(vars(Group,Repeat), ncol=11,nrow=2,strip.position = "bottom", scales = "free_x",shrink=TRUE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
		 geom_text(aes(label = round(r2,3)), 
            position = position_dodge(width = 0.9), size =3,vjust = -0.25)+
					 geom_text(aes(label = N,y=-0.005), 
            position = position_dodge(width = 0.9),size =2.5, vjust = -0.25)
			
            dev.off()

###############
#library(plyr)
library(ggplot2)
#library(reshape2)
library(grid)
a=read.csv('/Users/weig/Desktop/result_by_groups.csv')
b=read.csv('/Users/weig/Desktop/result_by_groups_sum.csv')
b$lower = b$r2_mean - b$r2_sd
b$upper = b$r2_mean + b$r2_sd
pd <- position_dodge(width = 0.2)
pdf(file =paste0('J:/PRS/PRS-CS_R2_by_sample_size.pdf'),width = 24,height = 18 )
par(mfrow=c(2,2))

for (data in c('unrel', 'ex.unrel')) {
  for (sex in c('women', 'men')) {

assign(paste0('gbase_',data,sex), ggplot(subset(b,Data==data &Sex==sex), aes(x=N_mean,y=r2_mean)) +
  geom_point(subset(a,Data==data &Sex==sex),position=pd,
             mapping = aes(x = N, y = r2),colour='grey') +  
  geom_point(subset(b,Data==data &Sex==sex),mapping = aes(x = N_mean, y = r2_mean),colour='red') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.3,colour='blue',linetype = "dashed") +
  scale_y_continuous(breaks = round(seq(min(a$r2), max(a$r2), by = 0.01),2)) +
  scale_x_continuous(breaks = round(seq(min(b$N_mean), max(b$N_mean), by = 1000),1)) +
  geom_line(subset(b,Data==data &Sex==sex),mapping = aes(x = N_mean, y = r2_mean),colour='black') +
  ggtitle(paste0('PRS-CS PRS R2 estimation of WHRI_RINT in ',gsub('ex.Training','Testing Data',gsub('unrel','Training Data',data)),' (',toupper(sex),')')) +
            theme(plot.title = element_text(hjust = 0.5)))

  }
}

pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(gbase_unrelwomen, vp = vplayout(1, 1))
print(gbase_unrelmen, vp = vplayout(1, 2))
print(gbase_ex.unrelwomen, vp = vplayout(2, 1))
print(gbase_ex.unrelmen, vp = vplayout(2, 2))

dev.off()


pdf(file =paste0('J:/PRS/PRS-CS_R2_by_sample_size.pdf'),width = 24,height = 18 )
par(mfrow=c(1,1))

ggplot(b, aes(x=N_mean,y=r2_mean)) +
  geom_point(a,position=pd,
             mapping = aes(x = N, y = r2),colour='grey') +  
  geom_point(b,mapping = aes(x = N_mean, y = r2_mean),colour='red') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.3,colour='blue',linetype = "dashed") +
  scale_y_continuous(breaks = round(seq(min(a$r2), max(a$r2), by = 0.01),2)) +
  scale_x_continuous(breaks = round(seq(min(b$N_mean), max(b$N_mean), by = 1000),1)) +
  geom_line(b,mapping = aes(x = N_mean, y = r2_mean),colour='black') +
  ggtitle('PRS-CS PRS R2 estimation of WHRI_RINT in Women') +
            theme(plot.title = element_text(hjust = 0.5))
dev.off()


# R2 1KG EAS CKB UKB LD references
result=as.data.frame(array(dim = c(1,8),NA))
for (ref in c('EAS','UKB','CKB')){
for (sub in c('unrel', 'ex.unrel')){
    for (sex in c('men','women')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',ref,'.',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in nrow(pheno)){
fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
result=rbind(result,c(ref,sub,n,sex,summary(fit)$r.squared,summary(fit)$coefficients[2,c(1,2,4)]))
}
}
}
}
result=result[-1,]    
names(result)=c('LD-Ref','Data','N','Sex','r2','BETA','SE','P')
result=result[,c('LD-Ref','Data','N','Sex','BETA','SE','P','r2')]
result

##-------------------
#   run PRS-CS using CKB ld references
for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/shared/ckb_bfile/ld/$chrom \
--extract /well/ckb/users/ahd580/tools/1KGPhase3_hapmap3_snp/w_hm3.snplist \
--make-bed \
--out /well/ckb/shared/ckb_bfile/ld/CKB_hm3$chrom
done

for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/shared/ckb_bfile/ld/CKB_hm3$chrom \
--freq \
--out /well/ckb/shared/ckb_bfile/ld/CKB_hm3_freq$chrom
done

ls /well/ckb/shared/ckb_bfile/ld/CKB_hm3*.bed | sed -e 's/\.bed//g' > /well/ckb/shared/ckb_bfile/ld/merge_list.txt

/well/ckb/users/ahd580/tools/plink \
--merge-list /well/ckb/shared/ckb_bfile/ld/merge_list.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW
library(data.table)
# Read in bim file for reference
bim_all<-NULL
for(i in 1:22){
    bim<-fread(paste0('/well/ckb/shared/ckb_bfile/ld/CKB_hm3',i,'.bim'))
    bim_all<-rbind(bim_all, bim)
}

# Read in the freq files to create the snpinfo file
frq_all<-NULL
for(i in 1:22){
    frq<-fread(paste0('/well/ckb/shared/ckb_bfile/ld/CKB_hm3_freq',i,'.frq'))
    frq_all<-rbind(frq_all, frq)
}

# Create snpinfo file
frq_bim<-merge(bim_all, frq_all[,c('SNP','MAF','NCHROBS')], by.x='V2', by.y='SNP')
names(frq_bim)<-c('SNP','CHR','POS','BP','A1','A2','MAF','NCHROBS')
frq_bim<-frq_bim[,c('CHR','SNP','BP','A1','A2','MAF','NCHROBS')]
frq_bim<-frq_bim[order(frq_bim$CHR, frq_bim$BP),]

# Remove SNPs with a MAF < 0.01
frq_bim<-frq_bim[frq_bim$MAF > 0.01,]
# Remove SNPs with a missingness > 0.01
frq_bim<-frq_bim[frq_bim$NCHROBS > (max(frq_bim$NCHROBS)*0.99),]
frq_bim$NCHROBS<-NULL

# Extract snplist for each range in LD Block bed files download https://bitbucket.org/nygcresearch/ldetect-data/src/master/ASN/fourier_ls-all.bed
bed<-fread('/well/ckb/users/ahd580/tools/1KGPhase3_hapmap3_snp/ASN_fourier_ls-all.bed')
bed$chr<-as.numeric(gsub('chr','',bed$chr))

blk_chr<-bed$chr
write.table(blk_chr,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/blk_chr'), col.names=F, row.names=F, quote=F)

blk_size<-NULL
snpinfo_1kg_hm3<-NULL
for(i in 1:dim(bed)[1]){
  frq_bim_i<-frq_bim[frq_bim$CHR == bed$chr[i] & frq_bim$BP > bed$start[i] & frq_bim$BP < bed$stop[i],]
  blk_size<-c(blk_size,dim(frq_bim_i)[1])
  snpinfo_1kg_hm3<-rbind(snpinfo_1kg_hm3,frq_bim_i)
  write.table(frq_bim_i$SNP,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/Block_',i,'.snplist'), col.names=F, row.names=F, quote=F)
}

write.table(blk_size,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/blk_size'), col.names=F, row.names=F, quote=F)

write.table(snpinfo_1kg_hm3, paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/snpinfo_1kg_hm3'), col.names=T, row.names=F, quote=F)

q()


# Then, we need to calculate LD within each block using PLINK.
   

#!/bin/bash
#files name PRSCS_ckbld_block_r2.sh
#$ -cwd -V -N CKBLD
#$ -q short.qc
#$ -P ckb.prjc
#$ -t 1-1445
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW \
--extract /well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/Block_${SGE_TASK_ID}.snplist \
--r square \
--out /well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW.Block_${SGE_TASK_ID}

qsub PRSCS_ckbld_block_r2.sh

for block in $(seq 1 1445);do 
if [ ! -f /well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW.Block_${block}.ld ]; then
echo $block
fi
done

# Block 557 doesn't contain any SNPs


R --vanilla

library(data.table)

for(i in c(1:556,558:1445)){
  LD_mat<-as.matrix(fread(paste0('/well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW.Block_',i,'.ld')))
  if(sum(!is.finite(LD_mat)) > 0){
    print(i)
  }
}

q()

export HDF5_USE_FILE_LOCKING='FALSE'
python /well/ckb/users/ahd580/PGS/zammy/script/write_ldblk_1KG_EUR.py

# run prscs

for i in $(seq 1 22);do
ln -s /well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/ldblk_ckb_chr${i}.hdf5 \
/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_CKB/LD_Blocks/ldblk_1kg_chr${i}.hdf5;
done 


qsub PRScs.CKB.sh men 203124

qsub PRScs.CKB.sh women 240145

for SEX in women men ;do
qsub PRScs.plinkscore.CKB.sh $SEX 
done


##-------------------
#   run PRS-CS using UKB 15K ld references
for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr$chrom \
--extract /well/ckb/users/ahd580/tools/1KGPhase3_hapmap3_snp/w_hm3.snplist \
--make-bed \
--out /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3$chrom
done

for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3$chrom \
--freq \
--out /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_freq$chrom
done

ls /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3*.bed | sed -e 's/\.bed//g' > /well/ckb/users/ahd580/tools/ukb15kref/merge_list.txt

/well/ckb/users/ahd580/tools/plink \
--merge-list /well/ckb/users/ahd580/tools/ukb15kref/merge_list.txt \
--make-bed \
--out /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW

R --vanilla 
library(data.table)
# Read in bim file for reference
bim_all<-NULL
for(i in 1:22){
    bim<-fread(paste0('/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3',i,'.bim'))
    bim_all<-rbind(bim_all, bim)
}

# Read in the freq files to create the snpinfo file
frq_all<-NULL
for(i in 1:22){
    frq<-fread(paste0('/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_freq',i,'.frq'))
    frq_all<-rbind(frq_all, frq)
}

# Create snpinfo file
frq_bim<-merge(bim_all, frq_all[,c('SNP','MAF','NCHROBS')], by.x='V2', by.y='SNP')
names(frq_bim)<-c('SNP','CHR','POS','BP','A1','A2','MAF','NCHROBS')
frq_bim<-frq_bim[,c('CHR','SNP','BP','A1','A2','MAF','NCHROBS')]
frq_bim<-frq_bim[order(frq_bim$CHR, frq_bim$BP),]

# Remove SNPs with a MAF < 0.01
frq_bim<-frq_bim[frq_bim$MAF > 0.01,]
# Remove SNPs with a missingness > 0.01
frq_bim<-frq_bim[frq_bim$NCHROBS > (max(frq_bim$NCHROBS)*0.99),]
frq_bim$NCHROBS<-NULL

# Extract snplist for each range in LD Block bed files
bed<-fread('/well/ckb/users/ahd580/tools/1KGPhase3_hapmap3_snp/EUR_fourier_ls-all.bed')
bed$chr<-as.numeric(gsub('chr','',bed$chr))

blk_chr<-bed$chr
write.table(blk_chr,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/blk_chr'), col.names=F, row.names=F, quote=F)

blk_size<-NULL
snpinfo_1kg_hm3<-NULL
for(i in 1:dim(bed)[1]){
  frq_bim_i<-frq_bim[frq_bim$CHR == bed$chr[i] & frq_bim$BP > bed$start[i] & frq_bim$BP < bed$stop[i],]
  blk_size<-c(blk_size,dim(frq_bim_i)[1])
  snpinfo_1kg_hm3<-rbind(snpinfo_1kg_hm3,frq_bim_i)
  write.table(frq_bim_i$SNP,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/Block_',i,'.snplist'), col.names=F, row.names=F, quote=F)
}

write.table(blk_size,paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/blk_size'), col.names=F, row.names=F, quote=F)

write.table(snpinfo_1kg_hm3, paste0('/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/snpinfo_1kg_hm3'), col.names=T, row.names=F, quote=F)

q()

# Then, we need to calculate LD within each block using PLINK.
   

#!/bin/bash
#files name PRSCS_ukbld_block_r2.sh
#$ -cwd -V -N CKBLD
#$ -q short.qc
#$ -P ckb.prjc
#$ -t 1-1703
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW \
--extract /well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/Block_${SGE_TASK_ID}.snplist \
--r square \
--out /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW.Block_${SGE_TASK_ID}

qsub  PRSCS_ukbld_block_r2.sh

for block in $(seq 1 1703);do 
if [ ! -f /well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW.Block_${block}.ld ]; then
echo $block
fi
done

# Block 90 doesn't contain any SNPs


R --vanilla

library(data.table)

for(i in c(1:89,91:1703)){
  LD_mat<-as.matrix(fread(paste0('/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW.Block_',i,'.ld')))
  if(sum(!is.finite(LD_mat)) > 0){
    print(i)
  }
}

q()

export HDF5_USE_FILE_LOCKING='FALSE'
python /well/ckb/users/ahd580/PGS/zammy/script/UKB_write_ldblk_1KG_EUR.py

for i in $(seq 1 22);do
ln -s /well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/ldblk_ukb_chr${i}.hdf5 \
/well/ckb/users/ahd580/tools/PRScs/PRScs_LD_matrix_UKB/LD_Blocks/ldblk_1kg_chr${i}.hdf5;
done 

qsub PRScs.UKB.sh men 203124

qsub PRScs.UKB.sh women 240145

for SEX in women men ;do
qsub PRScs.plinkscore.UKB.sh $SEX 
done

##################
#updated the /well/ckb/users/ahd580/tools/PRScs to the latest version 15-Oct-2020 which could deal with the non-atgc snps in the summary stats
#ckb ukb ld references using MAF 0.001 cutoff

/well/ckb/users/ahd580/PGS/zammy/script/PRScs.maf001_CKB_script.sh
/well/ckb/users/ahd580/PGS/zammy/script/PRScs.maf001_UKB_script.sh

awk 'NR==1||(length($2)==1&&length($3)==1)'  /well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.prscs.txt >/well/ckb/users/ahd580/PGS/zammy/WHR_men_SNP.bgen.prscs.snp.only.txt
awk 'NR==1||(length($2)==1&&length($3)==1)'  /well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen.prscs.txt >/well/ckb/users/ahd580/PGS/zammy/WHR_women_SNP.bgen.prscs.snp.only.txt

qsub PRScs.maf001_UKB.sh men 203124

qsub PRScs.maf001_UKB.sh women 240145

qsub PRScs.maf001_CKB.sh men 203124

qsub PRScs.maf001_CKB.sh women 240145

for SEX in women men ;do
qsub PRScs.maf001_plinkscore.CKB.sh $SEX 
done

for SEX in women men ;do
qsub PRScs.maf001_plinkscore.UKB.sh $SEX 
done

##R2 calculation
result=as.data.frame(array(dim = c(1,8),NA))
for (ref in c('UKB','CKB')){
for (sub in c('unrel', 'ex.unrel')){
    for (sex in c('men','women')){
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
for (i in 1:22 ) {
    score=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/maf001',ref,'.',sex,'/chr',i,'.score.profile'),header=T)
    names(score)[6]=paste0('chr',i)
    score=score[,c(1:2,6)]
    pheno=merge(pheno,score,by=c('FID','IID'),all.x=T)
}
pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
for (n in nrow(pheno)){
fit=lm(WHR_RINT~prs,data=pheno[ sample(nrow(pheno),n),])
result=rbind(result,c(ref,sub,n,sex,summary(fit)$r.squared,summary(fit)$coefficients[2,c(1,2,4)]))
}
}
}
}
result=result[-1,]    
names(result)=c('LD-Ref','Data','N','Sex','r2','BETA','SE','P')
result=result[,c('LD-Ref','Data','N','Sex','BETA','SE','P','r2')]
result


##Lassosum
for subdata in 'ex.unrel' 'unrel';do
for sex in women men;do
for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink --bfile /well/ckb/shared/ckb_bfile/chr${chrom} \
--keep /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.sample.txt \
--make-bed \
--out /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.chr${chrom}
done
done
done

for subdata in 'ex.unrel' 'unrel';do
for sex in women men;do
ls /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.chr*.bed | sed -e 's/\.bed//g' > /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.chrall.merge_list.txt

/well/ckb/users/ahd580/tools/plink \
--merge-list /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.chrall.merge_list.txt \
--make-bed \
--out /well/ckb/users/ahd580/PGS/zammy/PRScs/$subdata.ckb_${sex}.GW

done
done


### prepare the fam file with residuals of the WHR_RINT adjusted form array, PCs
for (sub in c('unrel', 'ex.unrel')){
        for (sex in c('women','men')){
            pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
            for (i in 1:22 ) {
                    fam=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/',sub,'.ckb_',sex,'.chr',i,'.fam'),header=F)
                        names(fam)=c('FID','IID','MID','PID','SEX','Pheno')
                        fam$order=1:(nrow(fam))
                            fam=merge(fam,pheno,by=c('FID','IID'),all.x=T)
                                fam$WHR_RINT_res=ifelse(is.na(fam$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=fam)))
                                    fam$Pheno=fam$WHR_RINT_res
                                    fam=fam[order(fam$order),]
                                        write.table(fam[,1:6],paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/',sub,'.ckb_',sex,'.chr',i,'.fam'),row.names=F,quote=F,col.names=F)
                                            }
                                                }
                                                    }


 cp /well/ckb/users/ahd580/PGS/zammy/lassosum/*.fam /well/ckb/users/ahd580/PGS/zammy/PRScs/

## limite snps to hapmap3 snps 
#update pheno colunm by chr and subset
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/ckb_all_pheno_gsid.txt'),header=T)
for (sub in c('unrel', 'ex.unrel')){
        for (sex in c('women','men')){
            pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
            for (i in 1:22 ) {
                    fam=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sub,'_',sex,'_chr',i,'.fam'),header=F)
                        names(fam)=c('FID','IID','MID','PID','SEX','Pheno')
                        fam$order=1:(nrow(fam))
                            fam=merge(fam,pheno,by=c('FID','IID'),all.x=T)
                                fam$WHR_RINT_res=ifelse(is.na(fam$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=fam)))
                                    fam$Pheno=fam$WHR_RINT_res
                                    fam=fam[order(fam$order),]
                                        write.table(fam[,1:6],paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sub,'_',sex,'_chr',i,'.fam'),row.names=F,quote=F,col.names=F)
                                            }
                                                }
                                                    }
													
													
#update pheno colunm by subset
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/ckb_all_pheno_gsid.txt'),header=T)
for (sub in c('unrel', 'ex.unrel')){
        for (sex in c('women','men')){
            pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
                    fam=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sub,'_',sex,'.fam'),header=F)
                        names(fam)=c('FID','IID','MID','PID','SEX','Pheno')
                        fam$order=1:(nrow(fam))
                            fam=merge(fam,pheno,by=c('FID','IID'),all.x=T)
                                fam$WHR_RINT_res=ifelse(is.na(fam$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=fam)))
                                    fam$Pheno=fam$WHR_RINT_res
                                    fam=fam[order(fam$order),]
                                        write.table(fam[,1:6],paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sub,'_',sex,'.fam'),row.names=F,quote=F,col.names=F)
                                            }
                                                }
												
													
#update pheno colunm by chr 													
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/ckb_all_pheno_gsid.txt'),header=T)
            for (i in 1:22 ) {
                    fam=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3',i,'.fam'),header=F)
                        names(fam)=c('FID','IID','MID','PID','SEX','Pheno')
                        fam$order=1:(nrow(fam))
                            fam=merge(fam,pheno,by=c('FID','IID'),all.x=T)
                                fam$WHR_RINT_res=ifelse(is.na(fam$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=fam)))
                                    fam$Pheno=fam$WHR_RINT_res
                                    fam=fam[order(fam$order),]
                                        write.table(fam[,1:6],paste0('/well/ckb/shared/ckb_bfile/CKB_hm3',i,'.fam'),row.names=F,quote=F,col.names=F)
										}
#update pheno colunm whole dataset			
pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/ckb_all_pheno_gsid.txt'),header=T)
fam=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW.fam'),header=F)
	names(fam)=c('FID','IID','MID','PID','SEX','Pheno')
	fam$order=1:(nrow(fam))
		fam=merge(fam,pheno,by=c('FID','IID'),all.x=T)
			fam$WHR_RINT_res=ifelse(is.na(fam$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=fam)))
				fam$Pheno=fam$WHR_RINT_res
				fam=fam[order(fam$order),]
					write.table(fam[,1:6],paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW.fam'),row.names=F,quote=F,col.names=F)

lassosum.R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sex=args[1]
N_sample=as.numeric(args[2])
ldfile=args[3]
pop=args[4]
ld=args[5]
library('lassosum')
setwd(system.file("data", package="lassosum"))
library(data.table)
### Read summary statistics file ###
ss <- fread(paste0("/well/ckb/users/ahd580/PGS/zammy/WHR_",sex,"_SNP.bgen.txt"))
head(ss)

### Specify the PLINK file stub of the reference panel ###
#ref.bfile <- paste0(ldfile,chrom)

### Specify the PLINK file stub of the test data ###
#test.bfile <- paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/','unrel.ckb_',sex,'.chr',chrom)

### Read LD region file ###
LDblocks <- pop # This will use LD regions as defined in Berisa and Pickrell (2015) for the Asian population and the hg19 genome.
# Other alternatives available. Type ?lassosum.pipeline for more details.

### SNP-wise correlations, Converted from p-values via the p2cor function ###
cor <- p2cor(p = ss$P_value, n = N_sample, sign=ss$beta)
library("fdrtool")
library(parallel)
cl <- makeCluster(6, type="FORK") # Parallel over 2 nodes
for (chrom in 1:22){
assign('ref.bfile',paste0(ldfile,chrom))
assign('test.bfile',paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/','unrel.ckb_',sex,'.chr',chrom))

assign(paste0('out',chrom),lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$BP,
                         A1=ss$a1, A2=ss$a2, # A2 is not required but advised
                         ref.bfile=ref.bfile, test.bfile=test.bfile,
                         LDblocks = LDblocks,cluster=cl))
}
out <- merge(out1, out2,out3, out4,out5, out6,out7, out8,out9, out10,out11, out12,out13, out14,out15, out16,out17, out18,out19, out20, out21, out22)

assign('test2.bfile',paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/','ex.unrel.ckb_',sex,'.GW'))
assign('test1.bfile',paste0('/well/ckb/users/ahd580/PGS/zammy/PRScs/','unrel.ckb_',sex,'.GW'))

v <- validate(out)
outnew <- subset(out, s=v$best.s, lambda=v$best.lambda)
v2 <- validate(outnew, test.bfile=test2.bfile)

vpseud<- pseudovalidate(out, test.bfile=test1.bfile)
outpseud <- subset(out, s=vpseud$best.s, lambda=vpseud$best.lambda)

vpseud1 <- validate(outpseud)
vpseud2 <- validate(outpseud, test.bfile=test2.bfile)

results=as.data.frame(array(NA,dim=c(4,8)))
results[1,]<- c('validation',pop,'unrel',sex,summary(lm(pheno~best.pgs,data=v$results.table))$coefficients[2,c(1,2,4)],summary(lm(pheno~best.pgs,data=v$results.table))$r.squared)
results[2,]<- c('validation',pop,'ex.unrel',sex,summary(lm(pheno~best.pgs,data=v2$results.table))$coefficients[2,c(1,2,4)],summary(lm(pheno~best.pgs,data=v2$results.table))$r.squared)
results[3,]<- c('Pseuo_validation',pop,'unrel',sex,summary(lm(pheno~best.pgs,data=vpseud1$results.table))$coefficients[2,c(1,2,4)],summary(lm(pheno~best.pgs,data=vpseud1$results.table))$r.squared)
results[4,]<- c('Pseuo_validation',pop,'ex.unrel',sex,summary(lm(pheno~best.pgs,data=vpseud2$results.table))$coefficients[2,c(1,2,4)],summary(lm(pheno~best.pgs,data=vpseud2$results.table))$r.squared)

write.table(results,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/',sex,ld,pop,'results'),quote=F,row.names=F)

write.table(out$sumstats,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','unrel.ckb_',sex,'.chrAll','snplist',ld,pop),quote=F,row.names=F)
write.table(v$best.beta,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','unrel.ckb_',sex,'.chrAll','beta',ld,pop),quote=F,row.names=F)
write.table(v$results.table,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','unrel.ckb_',sex,'.chrAll','table',ld,pop),quote=F,row.names=F)
write.table(v2$results.table,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','ex.unrel.ckb_',sex,'.chrAll','table',ld,pop),quote=F,row.names=F)
write.table(vpseud1$best.beta,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','unrel.ckb_',sex,'.chrAll','beta_pseud',ld,pop),quote=F,row.names=F)
write.table(vpseud1$results.table,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','unrel.ckb_',sex,'.chrAll','table_pseud',ld,pop),quote=F,row.names=F)
write.table(vpseud2$results.table,paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/','ex.unrel.ckb_',sex,'.chrAll','table_pseud',ld,pop),quote=F,row.names=F)


##
qsub lassosum.sh men 203124 /well/ckb/shared/ckb_bfile/ld/ ASN.hg19 CKB10k
qsub lassosum.sh women 240145 /well/ckb/shared/ckb_bfile/ld/ ASN.hg19 CKB10k

qsub lassosum.sh men 203124 /well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr EUR.hg19 UKB15k
qsub lassosum.sh women 240145 /well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr EUR.hg19 UKB15k


qsub lassosum_hm3.sh men 203124 /well/ckb/shared/ckb_bfile/ld/ ASN.hg19 CKB10k
qsub lassosum_hm3.sh women 240145 /well/ckb/shared/ckb_bfile/ld/ ASN.hg19 CKB10k

qsub lassosum_hm3.sh men 203124 /well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr EUR.hg19 UKB15k
qsub lassosum_hm3.sh women 240145 /well/ckb/users/ahd580/tools/ukb15kref/random15k_common.chr EUR.hg19 UKB15k

for (ref in c('UKB15kEUR.hg19','CKB10kASN.hg19')){
for (sub in c('unrel', 'ex.unrel')){
        for (sex in c('men','women')){
            pheno=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/',sub,'.ckb_',sex,'_pheno_gsid.txt'),header=T)
            for (chrom in 1:22 ) {
                    fam=read.table(paste0('/well/ckb/users/ahd580/PGS/zammy/lassosum/',sub,'.ckb_',sex,'.chr',chrom,'table',ref),header=T)
                    names(fam)[4]=paste0('best.pgs',chrom)
                    fam=fam[,c(1:2,4)]
                    pheno=merge(pheno,fam,by=c('FID','IID'),all.y=T)
            }
            pheno$prs=rowSums(pheno[,24:45])
pheno$WHR_RINT_res=ifelse(is.na(pheno$array),NA,residuals(lm(WHR_RINT~array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno)))
fit=lm(WHR_RINT~prs,data=pheno)
fit=lm(WHR_RINT~prs,data=pheno)
print (paste0(ref,'_',sub,'_',sex,' ',summary(fit)$r.squared))
print (summary(fit)$coefficients[2,c(1,2,4)])
}
}
}


##ldpred2.R
qrsh -q short.qc -P ckb.prjc -pe shmem 4

##1 merge chromose bfiles 


for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/shared/ckb_bfile/chr$chrom \
--extract /well/ckb/users/ahd580/tools/1KGPhase3_hapmap3_snp/w_hm3.snplist \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3$chrom &
done

for chrom in $(seq 1 22);do 
/well/ckb/users/ahd580/tools/plink --bfile \
/well/ckb/shared/ckb_bfile/CKB_hm3$chrom \
--freq \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_freq$chrom
done

ls /well/ckb/shared/ckb_bfile/CKB_hm3*.bed | sed -e 's/\.bed//g' > /well/ckb/shared/ckb_bfile/merge_list.txt

/well/ckb/users/ahd580/tools/plink \
--merge-list /well/ckb/shared/ckb_bfile/merge_list.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_women_pheno_gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_women

ls -r /well/ckb/shared/ckb_bfile/CKB_hm3_GW_*unrel_men.fam|xargs cat|awk '$6!="NA" {print $1,$2}' > /well/ckb/shared/ckb_bfile/CKB_hm3_GW_men.gsid.txt
ls -r /well/ckb/shared/ckb_bfile/CKB_hm3_GW_*unrel_women.fam|xargs cat|awk '$6!="NA" {print $1,$2}' > /well/ckb/shared/ckb_bfile/CKB_hm3_GW_women.gsid.txt

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/shared/ckb_bfile/CKB_hm3_GW_men.gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_men

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/shared/ckb_bfile/CKB_hm3_GW_women.gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_women

for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_women \
--chr $chrom \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_women_chr$chrom
done

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/users/ahd580/PGS/zammy/unrel.ckb_women_pheno_gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_women

for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_women \
--chr $chrom \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_women_chr$chrom
done

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/users/ahd580/PGS/zammy/ex.unrel.ckb_men_pheno_gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_men

for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_men \
--chr $chrom \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_ex.unrel_men_chr$chrom
done

/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW \
--keep /well/ckb/users/ahd580/PGS/zammy/unrel.ckb_men_pheno_gsid.txt \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_men

for chrom in $(seq 1 22);do
/well/ckb/users/ahd580/tools/plink \
--bfile /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_men \
--chr $chrom \
--make-bed \
--out /well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_men_chr$chrom
done

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sex=args[1]
N_sample=as.numeric(args[2])
NCORES <-as.numeric(args[3])


library(bigsnpr)
library(bigreadr)
library(data.table)
library(bigsparser)


a=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sex,'.fam'),header=F)
b=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_',sex,'.fam'),header=F)
c=merge(a,b,by='V1',all.x=T)
length(row.names(c[!is.na(c$V6.y),]))
ind.val<-as.numeric(row.names(c[!is.na(c$V6.y),]))
ind.test<-as.numeric(row.names(c[is.na(c$V6.y),]))
#target genotype file 
test <- snp_attach(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sex,'.rds'))
names(test$map)= c("chr", "rsid", "genetic.dist", "pos", "a0", "a1")
map_test<-test$map[-(2:3)]
names(map_test) <- c("chr", "pos", "a0", "a1")

## Information for the variants provided in the LD reference
map_ldref <- readRDS("/well/ckb/users/ahd580/tools/ldpred2/map.rds")
## summary statistics
sumstats <- fread2(paste0("/well/ckb/users/ahd580/PGS/zammy/WHR_",sex,"_SNP.bgen.txt"))
sumstats <- sumstats[,c("CHR", "BP", "a1", "a2","beta","se",'P_value')]
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se",'P')
sumstats$n_eff <- N_sample
#match snp a1 a0 with summary stats, target map and ld ref
df_beta <- snp_match(sumstats, map_ldref)
(df_beta <- tidyr::drop_na(tibble::as_tibble(df_beta)))
df_beta2<-snp_match(sumstats, map_test)
(df_beta2 <- tidyr::drop_na(tibble::as_tibble(df_beta2)))
# Here, you also want to restrict to the variants present
# in your test data as well. For this, you can use something like

in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], df_beta2[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]
in_test <- vctrs::vec_in(df_beta2[, c("chr", "pos")], df_beta[, c("chr", "pos")])
df_beta2 <- df_beta2[in_test, ]
dir<-cbind(df_beta,df_beta2)
dir$dir=ifelse(dir[,5]==dir[,21],1,-1)

tmp <- tempfile(tmpdir = "tmp-data")

#load UKB LD matrices

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("/well/ckb/users/ahd580/tools/ldpred2/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
       ld <- Matrix::colSums(corr_chr^2)
    corr <- as_SFBM(corr_chr, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr_chr^2))
    corr$add_columns(corr_chr, nrow(corr))
  }
}


# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(df_beta),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))

h2_est <- ldsc[["h2"]]

G <- test$genotypes
y=test$fam

#infinitesimal validation data
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
beta_inf_target<-beta_inf*dir$dir

                    
# LDpred2-grid
  (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  beta_grid_target<-beta_grid*dir$dir
  params$sparsity <- colMeans(beta_grid_target == 0)

  bigparallelr::set_blas_ncores(NCORES)
  pred_grid <- big_prodMat(G, beta_grid_target, ind.row = ind.val ,
                           ind.col = df_beta2[["_NUM_ID_"]])

  params$score <- big_univLinReg(as_FBM(pred_grid), y[ind.val,'affection'])$score

  library(dplyr)
  best_beta_grid_nosp <- params %>%
    mutate(id = row_number()) %>%
    filter(!sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]

  best_beta_grid_sp <- params %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
	
	# LDpred2-auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.9, 30),
                                 ncores = NCORES)
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  beta_auto_target<-beta_auto*dir$dir

  #auto validation data
  pred_auto_val <- big_prodMat(G, beta_auto_target, ind.row = ind.val ,
                           ind.col = df_beta2[["_NUM_ID_"]])
  sc <- apply(pred_auto_val, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto_val <- rowMeans(beta_auto_target[, keep])

  #auto test data
  pred_auto_test <- big_prodMat(G, beta_auto_target, ind.row = ind.test ,
                           ind.col = df_beta2[["_NUM_ID_"]])
  sc <- apply(pred_auto_test, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto_test <- rowMeans(beta_auto_target[, keep])



  # compute predictions for val set
  betas_val <- cbind(beta_inf_target, final_beta_auto_val)

  pred_val <- big_prodMat(G, betas_val, ind.row = ind.val,
  		   ind.col = df_beta2[["_NUM_ID_"]])
result<-NULL
lm_val<-cbind(y[ind.val,],pred_val)
names(lm_val)[7:8]=c('beta_inf_target', 'final_beta_auto_val')
for (i in 7:8) {
glm<-summary(lm(lm_val[,'affection']~lm_val[,i]))
result<-rbind(result,c('val',names(lm_val)[i],glm$coefficients[2,c(1,2,4)],glm$r.squared))
}

  # compute predictions for test set
  betas_test <- cbind(beta_inf_target, best_beta_grid_nosp, best_beta_grid_sp,
                 final_beta_auto_test)
  pred_test <- big_prodMat(G, betas_test, ind.row = ind.test,
  		   ind.col = df_beta2[["_NUM_ID_"]])

lm_test<-cbind(y[ind.test,],pred_test)
names(lm_test)[7:10]=c('beta_inf_target', 'best_beta_grid_nosp', 'best_beta_grid_sp','final_beta_auto_test')
for (i in 7:10) {
glm<-summary(lm(lm_test[,'affection']~lm_test[,i]))
result<-rbind(result,c('test',names(lm_test)[i],glm$coefficients[2,c(1,2,4)],glm$r.squared))
}


  # save results
  res_val <- list(pred = setNames(as.data.frame(pred_val), colnames(betas_val)),
              params = params, auto = multi_auto[keep])
  saveRDS(res_val, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.res_val_",sex))


  res_test <- list(pred = setNames(as.data.frame(pred_test), colnames(betas_test)),
              params = params, auto = multi_auto[keep])
  saveRDS(res_test, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.res_test",sex))

  write.table(betas_val,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.betas_val",sex),row.names=F,quote=F)
  write.table(betas_test,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.betas_test",sex),row.names=F,quote=F)

  write.table(pred_val,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.pred_val",sex),row.names=F,quote=F)
  write.table(pred_test,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.pred_test",sex),row.names=F,quote=F)

  write.table(result,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/provided_ld.pred_test",sex,'result.txt'),row.names=F,quote=F)


#!/bin/bash
#ldpred2_provided_ld.sh
#$ -cwd -V -N ldpred2_provided_ld
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 12
gender=$1
N=$2
/well/ckb-share/local/bin/Rscript /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_provided_ld.R $gender $N ${NSLOTS:-1}

qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_provided_ld.sh women 240145
qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_provided_ld.sh men 203124


#create ld matrices uisng ukb and ckb data

args = commandArgs(trailingOnly=TRUE)

NCORES <-as.numeric(args[1])

library(bigsnpr)
library(bigreadr)
library(data.table)
library(bigsparser)
ukbld <- snp_attach("/well/ckb/users/ahd580/tools/ukb15kref/UKB_hm3_GW.rds")
G <- ukbld$genotypes
CHR <- as.integer(ukbld$map$chromosome)
POS <- ukbld$map$physical.pos

ind.val<-1:nrow(G)

#down load  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
#cd /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/
#for i in $(seq 1 22); do mv CEU-${i}-final.genetic.map chr$i.OMNI.interpolated_genetic_map;done

POS2 <- snp_asGeneticPos(CHR, POS, dir = "/well/ckb/users/ahd580/tools/gctb_2.02_Linux/CEU/", ncores = NCORES)

bigassertr::assert_dir("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ukb_corr")

for (chr in 1:22) {

  ind.chr <- which(CHR == chr)

  runonce::save_run(
    snp_cor(
      G, ind.row = ind.val, ind.col = ind.chr,
      alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
      ncores = NCORES
    ),
    file = paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ukb_corr/chr", chr, ".rds")
  )
}
  saveRDS(ukbld$map, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ukb_corr/map.rds"))



args = commandArgs(trailingOnly=TRUE)

NCORES <-as.numeric(args[1])

library(bigsnpr)
library(bigreadr)
library(data.table)
library(bigsparser)

ckbld <- snp_attach("/well/ckb/shared/ckb_bfile/ld/CKB_hm3_GW.rds")
G <- ckbld$genotypes
CHR <- as.integer(ckbld$map$chromosome)
POS <- ckbld$map$physical.pos

ind.val<-1:nrow(G)

# cd /well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/
# for i in $(seq 1 22); do mv CHB-${i}-final.genetic.map chr$i.OMNI.interpolated_genetic_map;done

POS2 <- snp_asGeneticPos(CHR, POS, dir = "/well/ckb/users/ahd580/tools/gctb_2.02_Linux/CHB/", ncores = NCORES)

bigassertr::assert_dir("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_corr")

for (chr in 1:22) {

  ind.chr <- which(CHR == chr)

  runonce::save_run(
    snp_cor(
      G, ind.row = ind.val, ind.col = ind.chr,
      alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
      ncores = NCORES
    ),
    file = paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_corr/chr", chr, ".rds")
  )
}

saveRDS(ckbld$map, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_corr/map.rds"))


#!/bin/bash
#ckb_ld_matrices.sh
#$ -cwd -V -N ckb_ld_matrices
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 12

/well/ckb-share/local/bin/Rscript ckb_ld_matrices.r ${NSLOTS:-1}


#!/bin/bash
# ukb_ld_matrices.sh
#$ -cwd -V -N ukb_ld_matrices
#$ -q short.qc
#$ -P ckb.prjc
#$ -pe shmem 12

/well/ckb-share/local/bin/Rscript ukb_ld_matrices.r

########
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sex=args[1]
N_sample=as.numeric(args[2])
NCORES <- as.numeric(args[3])

library(bigsnpr)
library(bigreadr)
library(data.table)
library(bigsparser)


a=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sex,'.fam'),header=F)
b=read.table(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_unrel_',sex,'.fam'),header=F)
c=merge(a,b,by='V1',all.x=T)
length(row.names(c[!is.na(c$V6.y),]))
ind.val<-as.numeric(row.names(c[!is.na(c$V6.y),]))
ind.test<-as.numeric(row.names(c[is.na(c$V6.y),]))
#target genotype file 
test <- snp_attach(paste0('/well/ckb/shared/ckb_bfile/CKB_hm3_GW_',sex,'.rds'))
names(test$map)= c("chr", "rsid", "genetic.dist", "pos", "a0", "a1")
map_test<-test$map[-(2:3)]
names(map_test) <- c("chr", "pos", "a0", "a1")

## Information for the variants provided in the LD reference
map_ldref <- readRDS("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_corr/map.rds")
map_ldref<-map_ldref[,-3]
names(map_ldref)= c("chr", "rsid", "pos","a1", "a0")
## summary statistics
sumstats <- fread2(paste0("/well/ckb/users/ahd580/PGS/zammy/WHR_",sex,"_SNP.bgen.txt"))
sumstats <- sumstats[,c("CHR", "BP", "a1", "a2","beta","se",'P_value')]
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se",'P')
sumstats$n_eff <- N_sample
#match snp a1 a0 with summary stats, target map and ld ref
df_beta <- snp_match(sumstats, map_ldref)
(df_beta <- tidyr::drop_na(tibble::as_tibble(df_beta)))
df_beta2<-snp_match(sumstats, map_test)
(df_beta2 <- tidyr::drop_na(tibble::as_tibble(df_beta2)))
# Here, you also want to restrict to the variants present
# in your test data as well. For this, you can use something like

in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], df_beta2[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]
in_test <- vctrs::vec_in(df_beta2[, c("chr", "pos")], df_beta[, c("chr", "pos")])
df_beta2 <- df_beta2[in_test, ]
dir<-cbind(df_beta,df_beta2)
dir$dir=ifelse(dir[,5]==dir[,21],1,-1)

tmp <- tempfile(tmpdir = "tmp-data")

#load CKB LD matrices

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
      ld <- Matrix::colSums(corr_chr^2)
    corr <- as_SFBM(corr_chr, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr_chr^2))
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(df_beta),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))

h2_est <- ldsc[["h2"]]

G <- test$genotypes
y=test$fam

#infinitesimal validation data
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
beta_inf_target<-beta_inf*dir$dir

                    
# LDpred2-grid
  (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  beta_grid_target<-beta_grid*dir$dir
  params$sparsity <- colMeans(beta_grid_target == 0)

  bigparallelr::set_blas_ncores(NCORES)
  pred_grid <- big_prodMat(G, beta_grid_target, ind.row = ind.val ,
                           ind.col = df_beta2[["_NUM_ID_"]])

  params$score <- big_univLinReg(as_FBM(pred_grid), y[ind.val,'affection'])$score

  library(dplyr)
  best_beta_grid_nosp <- params %>%
    mutate(id = row_number()) %>%
    filter(!sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]

  best_beta_grid_sp <- params %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
	
	# LDpred2-auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.9, 30),
                                 ncores = NCORES)
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  beta_auto_target<-beta_auto*dir$dir

  #auto validation data
  pred_auto_val <- big_prodMat(G, beta_auto_target, ind.row = ind.val ,
                           ind.col = df_beta2[["_NUM_ID_"]])
  sc <- apply(pred_auto_val, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto_val <- rowMeans(beta_auto_target[, keep])

  #auto test data
  pred_auto_test <- big_prodMat(G, beta_auto_target, ind.row = ind.test ,
                           ind.col = df_beta2[["_NUM_ID_"]])
  sc <- apply(pred_auto_test, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto_test <- rowMeans(beta_auto_target[, keep])



  # compute predictions for val set
  betas_val <- cbind(beta_inf_target, final_beta_auto_val)

  pred_val <- big_prodMat(G, betas_val, ind.row = ind.val,
  		   ind.col = df_beta2[["_NUM_ID_"]])
result<-NULL
lm_val<-cbind(y[ind.val,],pred_val)
names(lm_val)[7:8]=c('beta_inf_target', 'final_beta_auto_val')
for (i in 7:8) {
glm<-summary(lm(lm_val[,'affection']~lm_val[,i]))
result<-rbind(result,c('val',names(lm_val)[i],glm$coefficients[2,c(1,2,4)],glm$r.squared))
}

  # compute predictions for test set
  betas_test <- cbind(beta_inf_target, best_beta_grid_nosp, best_beta_grid_sp,
                 final_beta_auto_test)
  pred_test <- big_prodMat(G, betas_test, ind.row = ind.test,
  		   ind.col = df_beta2[["_NUM_ID_"]])
result<-NULL
lm_test<-cbind(y[ind.test,],pred_test)
names(lm_test)[7:10]=c('beta_inf_target', 'best_beta_grid_nosp', 'best_beta_grid_sp','final_beta_auto_test')
for (i in 7:10) {
glm<-summary(lm(lm_test[,'affection']~lm_test[,i]))
result<-rbind(result,c('test',names(lm_test)[i],glm$coefficients[2,c(1,2,4)],glm$r.squared))
}


  # save results
  res_val <- list(pred = setNames(as.data.frame(pred_val), colnames(betas_val)),
              params = params, auto = multi_auto[keep])
  saveRDS(res_val, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.res_val_",sex))


  res_test <- list(pred = setNames(as.data.frame(pred_test), colnames(betas_test)),
              params = params, auto = multi_auto[keep])
  saveRDS(res_test, paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.res_test",sex))

  write.table(betas_val,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.betas_val",sex),row.names=F,quote=F)
  write.table(betas_test,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.betas_test",sex),row.names=F,quote=F)

  write.table(pred_val,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.pred_val",sex),row.names=F,quote=F)
  write.table(pred_test,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.pred_test",sex),row.names=F,quote=F)

  write.table(result,paste0("/well/ckb/users/ahd580/PGS/zammy/ldpred2/ckb_ld.pred_test",sex,'result.txt'),row.names=F,quote=F)

qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_ckb_ld.sh women 240145
qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_ckb_ld.sh men 203124

qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_ukb_ld.sh women 240145
qsub /well/ckb/users/ahd580/PGS/zammy/script/ldpred2_ukb_ld.sh men 203124

