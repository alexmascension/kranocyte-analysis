BASEDIR="/mnt/data/Proyectos/kranocito/data/ARAUZO_03"

# Download CB Whitelist 
wget https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz -P ${BASEDIR}/CB_whitelist
gzip -d ${BASEDIR}/CB_whitelist/3M-february-2018.txt.gz

# Download mm39 data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz -P ${BASEDIR}/database

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz -P ${BASEDIR}/database

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz -P ${BASEDIR}/database

gzip -d ${BASEDIR}/database/*.gz

# Determine the sjdbOverhang parameter (Its L - 1)
zcat ${BASEDIR}/20230623/FASTQ_renamed/AZ7845_L001_R2_001.fastq.gz | awk 'NR%4==2 && ++n<=10000{len[length($0)]++} END{for (l in len) print l, len[l]}' | sort -n

# Create STAR index
docker pull gcfntnu/star:2.7.11b

docker run --rm -u $(id -u):$(id -g) -v ${BASEDIR}/database:/database -w /database \
  gcfntnu/star:2.7.11b \
  STAR --runThreadN 16 --runMode genomeGenerate \
  --genomeDir /database/mm39_STAR \
  --genomeFastaFiles /database/GCF_000001635.27_GRCm39_genomic.fna \
  --sjdbGTFfile /database/GCF_000001635.27_GRCm39_genomic.gtf \
  --sjdbOverhang 89 


# Run STARsolo

mkdir -p "$PWD/results_STAR"

for name in AZ7845 AZ7846; do
  echo "Running file with name ${name}."
  mkdir -p "$PWD/results_STAR/${name}"
  
  docker run --rm -u "$(id -u)":"$(id -g)" \
    -v "${BASEDIR}/database:/database" \
    -v "${BASEDIR}/20230623/FASTQ_renamed:/fastq" \
    -v "$PWD/CB_whitelist:/whitelists" \
    -v "$PWD/results_STAR:/results" \
    -w /results \
    gcfntnu/star:2.7.11b \
    STAR \
      --runThreadN 16 \
      --genomeDir /database/mm39_STAR \
      --readFilesIn \
        "/fastq/${name}_L001_R2_001.fastq.gz,/fastq/${name}_L002_R2_001.fastq.gz,/fastq/${name}_L003_R2_001.fastq.gz,/fastq/${name}_L004_R2_001.fastq.gz" \
        "/fastq/${name}_L001_R1_001.fastq.gz,/fastq/${name}_L002_R1_001.fastq.gz,/fastq/${name}_L003_R1_001.fastq.gz,/fastq/${name}_L004_R1_001.fastq.gz" \
      --readFilesCommand zcat \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist /whitelists/3M-february-2018.txt \
      --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
      --clipAdapterType CellRanger4 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloCellFilter EmptyDrops_CR \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --soloMultiMappers EM \
      --outFileNamePrefix "/results/${name}/"
done
