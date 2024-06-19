# This process is the downstream analysis of the MAGs generated from assembly and binning

# 1. strain-level MAGs dereplication per sample

for b in hostsname; do
mkdir /mnt/z/xb/metagenome/finalbins/${b}/straindRep
cd /mnt/z/xb/metagenome/finalbins/${b}/bins
ls *.fa > bins.list
dRep dereplicate \
	-p 64 \
	-g bins.list \
	-pa 0.9 \
	-sa 0.99 \
	-nc 0.30 \
	-comp 80 \
	-con 5 /mnt/z/xb/metagenome/finalbins/${b}/80straindRep && echo "** Strain level Dereplication done **"
done

# 2. species-level MAGs dereplication per sample (SGBs)

for b in hostsname; do
mkdir /mnt/z/xb/metagenome/finalbins/${b}/SGBdRep
cd /mnt/z/xb/metagenome/finalbins/${b}/straindRep/dereplicated_genomes
ls *.fa > bins.list
dRep dereplicate \
	-p 64 \
	-g bins.list \
	-pa 0.9 \
	-sa 0.95 \
	-nc 0.30 \
	-comp 80 \
	-con 5 /mnt/z/xb/metagenome/finalbins/${b}/SGBdRep && echo "** SGB cluster done **"
done

# 3. species-level MAGs dereplication in all samples (SGBs)

mkdir /mnt/z/xb/metagenome/finalbins/allsamples
for b in hostsname; do
cp -r /mnt/z/xb/metagenome/finalbins/${b}/straindRep/dereplicated_genomes/*.fa /mnt/z/xb/metagenome/finalbins/allsamples
done
cd /mnt/z/xb/metagenome/finalbins/allsamples
mkdir /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep
ls *.fa > bins.list
dRep dereplicate \
	-p 64 \
	-g bins.list \
	-pa 0.9 \
	-sa 0.95 \
	-nc 0.30 \
	-comp 80 \
	-con 5 /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep && echo "** SGB cluster done **"



