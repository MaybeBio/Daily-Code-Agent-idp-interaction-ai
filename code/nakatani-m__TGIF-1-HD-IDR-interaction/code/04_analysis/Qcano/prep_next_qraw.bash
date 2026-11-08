#!/bin/bash
#$ -l h_rt=10:00:00
#$ -S /bin/bash
#$ -l cpu_8=1
#$ -cwd

path=`pwd`
echo $path
STAGE=${path##*/}
STAGE_fin=65

weight_factor=1.0
PREV_STAGE=$(( $STAGE -1 ))
echo ${PREV_STAGE}
export WORKDIR=/gs/bs/hp170020/tkhs
export KKTOOLS=${WORKDIR}/local/kktools
OMEGAROOT=/gs/bs/hp170020/tkhs/local/omegagene/og0_v059n
OMEGATK=${OMEGAROOT}/toolkit

rm params_next.dat

FN_QRAWIS="qrawis_list.txt"
FN_WEIGHT="qrawis_weight.txt"
FN_QRAW="qraw_list.txt"
test_weight=0.95
weight=1.0
rm ${FN_QRAWIS}
rm ${FN_QRAW}
rm ${FN_WEIGHT}

for st in `seq ${STAGE} -1 1`
    do
    nser=`ls -1 ../${st}/n*/qraw_is.dat | wc -l`
    ls -1 ../${st}/n*/qraw_is.dat >> ${FN_QRAWIS}
    ls -1 ../${st}/n*/qraw.dat >> ${FN_QRAW}
    if [ $st == $STAGE_fin ]
    then
          for se in `seq 1 $nser`
          do
                  awk "BEGIN { print $weight }"
                  awk "BEGIN { print $weight }" >> ${FN_WEIGHT}
          done
          weight=`echo "scale=5; ${weight} * ${weight_factor}" | bc`
    else
          for se in `seq 1 $nser`
          do
                  awk "BEGIN { print $test_weight**($STAGE_fin - $st) }"
                  awk "BEGIN { print $test_weight**($STAGE_fin - $st) }" >> ${FN_WEIGHT}
          done
          weight=`echo "scale=5; ${test_weight} * ${weight_factor}" | bc`
    fi
done

python ${OMEGATK}/vcmd_sum_qraw_is.py --i_qraw_list ${FN_QRAWIS}  --out qraw_sum.dat  --i_weight ${FN_WEIGHT}
${OMEGAROOT}/qcano_szfit --i-cfg qcano_szfit.txt

python ${OMEGATK}/merge_vcmd_qraw.py -i qw_opt.dat -o params_next.dat --p-count 0 --i-qraw-list ${FN_QRAW} --o-qraw vcmd_next_qraw.dat  --i-weight ${FN_WEIGHT}

source /gs/bs/hp170020/nak/virtual_env/venv2507/bin/activate
python ../weighted_ensemble/power_params.py
deactivate
