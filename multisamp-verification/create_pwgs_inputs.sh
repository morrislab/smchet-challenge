#DO51958 PRAD-UK SP114941        0bfd1043-8170-e3e4-e050-11ac0c4860c5
#DO51958 PRAD-UK SP114937        0bfd1043-8172-e3e4-e050-11ac0c4860c5
#DO51958 PRAD-UK SP114944        0bfd1068-3fc3-a95b-e050-11ac0c4860c3
#DO51958 PRAD-UK SP114943        0bfd1043-8175-e3e4-e050-11ac0c4860c5
#DO51958 PRAD-UK SP114945        0bfd1043-8173-e3e4-e050-11ac0c4860c5

#DO51965 PRAD-UK SP114991        0bfd1068-3fdf-a95b-e050-11ac0c4860c3
#DO51965 PRAD-UK SP114992        0bfd1068-3fe1-a95b-e050-11ac0c4860c3
#DO51965 PRAD-UK SP114994        0bfebf9f-c77d-e57d-e050-11ac0d487827
#DO51965 PRAD-UK SP114987        0bfd1043-8181-e3e4-e050-11ac0c4860c5
#DO51965 PRAD-UK SP114996        0bfd1043-8180-e3e4-e050-11ac0c4860c5


set -euo pipefail

BASE_DIR=~/work/exultant-pistachio
RUN_CNV_PATH=$BASE_DIR/data/cnvs/runf
VCF_PATH=$BASE_DIR/data/variants/consensus
SAMPLE_SIZE=5000
SEX=male
PARSER_PATH=~/.apps/phylowgs/parser

PATIENT=patient2

#S1=0bfd1043-8170-e3e4-e050-11ac0c4860c5
#S2=0bfd1043-8172-e3e4-e050-11ac0c4860c5
#S3=0bfd1068-3fc3-a95b-e050-11ac0c4860c3
#S4=0bfd1043-8175-e3e4-e050-11ac0c4860c5
#S5=0bfd1043-8173-e3e4-e050-11ac0c4860c5

S1=0bfd1068-3fdf-a95b-e050-11ac0c4860c3
S2=0bfd1068-3fe1-a95b-e050-11ac0c4860c3
S3=0bfebf9f-c77d-e57d-e050-11ac0d487827
S4=0bfd1043-8181-e3e4-e050-11ac0c4860c5
S5=0bfd1043-8180-e3e4-e050-11ac0c4860c5

echo python2 $PARSER_PATH/create_phylowgs_inputs.py \
    --cnvs S1=$RUN_CNV_PATH/$S1.txt \
    --cnvs S2=$RUN_CNV_PATH/$S2.txt \
    --cnvs S3=$RUN_CNV_PATH/$S3.txt \
    --cnvs S4=$RUN_CNV_PATH/$S4.txt \
    --cnvs S5=$RUN_CNV_PATH/$S5.txt \
    --output-variants ${PATIENT}.sampled.ssm \
    --output-cnvs ${PATIENT}.sampled.cnv \
    --nonsubsampled-variants ${PATIENT}.nonsampled.ssm \
    --nonsubsampled-variants-cnvs ${PATIENT}.nonsampled.cnv \
    --sample-size $SAMPLE_SIZE \
    --vcf-type S1=pcawg_consensus \
    --vcf-type S2=pcawg_consensus \
    --vcf-type S3=pcawg_consensus \
    --vcf-type S4=pcawg_consensus \
    --vcf-type S5=pcawg_consensus \
    --sex $SEX \
    --verbose \
    S1=$VCF_PATH/${S1}.consensus.20160830.somatic.snv_mnv.vcf.gz \
    S2=$VCF_PATH/${S2}.consensus.20160830.somatic.snv_mnv.vcf.gz \
    S3=$VCF_PATH/${S3}.consensus.20160830.somatic.snv_mnv.vcf.gz \
    S4=$VCF_PATH/${S4}.consensus.20160830.somatic.snv_mnv.vcf.gz \
    S5=$VCF_PATH/${S5}.consensus.20160830.somatic.snv_mnv.vcf.gz \
    ">$PATIENT.stdout 2>$PATIENT.stderr"
