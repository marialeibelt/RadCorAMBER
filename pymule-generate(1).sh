#!/bin/bash
build_piece () {
    for i in $@ ; do
        echo input/$i
    done
}
generate() {
    order="$1"
    flag="$2"
    xi="$3"
    seed="$4"
    tgteff=1e-3
    dmax=0.0045

    if [[ "$order" == 0 ]] ; then
        npieces=1
        pieces=$(build_piece mp2mp0_e-p-_S0000048034X1.00000D1.00000_ITMX020x010.0M.mcmule)
        nenter=5000
        itmx=50
        tgteff=-1
    elif [[ "$order" == 1 ]] ; then
        npieces=3
        nenter=5000
        itmx=50
        nsub=30

        if [[ "$xi" == "01" ]] ; then
            pieces=$(build_piece \
                mp2mpNLO0_e-p-_S0000044858X0.10000D0.10000_ITMX020x010.0M.mcmule \
                mp2mpR15_e-p-_S0000011656X0.10000D0.10000_ITMX020x050.0M.mcmule \
                mp2mpR35_e-p-_S0000054321X0.10000D0.10000_ITMX020x050.0M.mcmule)
        elif [[ "$xi" == "10" ]] ; then
            pieces=$(build_piece \
                mp2mpNLO0_e-p-_S0000072707X1.00000D1.00000_ITMX020x010.0M.mcmule \
                mp2mpR15_e-p-_S0000016843X1.00000D1.00000_ITMX020x050.0M.mcmule \
                mp2mpR35_e-p-_S0000050947X1.00000D1.00000_ITMX020x050.0M.mcmule)
        fi
    elif [[ "$order" == 2 ]] ; then
        npieces=5
        nenter=40000
        itmx=80

        if [[ "$xi" == "01" ]] ; then
            pieces=$(build_piece \
                mp2mpNNLO0_e-p-_S0001194496X0.10000D0.10000_ITMX020x010.0M.mcmule \
                mp2mpNNLO115_e-p-_S0001155921X0.10000D0.10000_ITMX020x050.0M.mcmule \
                mp2mpNNLO135_e-p-_S0001139033X0.10000D0.10000_ITMX020x050.0M.mcmule \
                mp2mpRR1516_e-p-_S0000079013X0.10000D0.10000_ITMX020x100.0M.mcmule \
                mp2mpRR3536_e-p-_S0000021942X0.10000D0.10000_ITMX020x100.0M.mcmule )
        elif [[ "$xi" == "10" ]] ; then
            pieces=$(build_piece \
                mp2mpNNLO0_e-p-_S0001165308X1.00000D1.00000_ITMX020x010.0M.mcmule \
                mp2mpNNLO115_e-p-_S0001112439X1.00000D1.00000_ITMX020x050.0M.mcmule \
                mp2mpNNLO135_e-p-_S0001141606X1.00000D1.00000_ITMX020x050.0M.mcmule \
                mp2mpRR1516_e-p-_S0000021884X1.00000D1.00000_ITMX020x100.0M.mcmule \
                mp2mpRR3536_e-p-_S0000064553X1.00000D1.00000_ITMX020x100.0M.mcmule )
        fi
    fi

    if [[ "$flag" == "N" ]] ; then
        features="V"
    elif [[ "$flag" == "C" ]] ; then
        features="VC\n$dmax"
    elif [[ "$flag" == "S" ]] ; then
        features="VCS\n$dmax\n$nsub"
    elif [[ "$flag" == "U" ]] ; then
        features="VU\n$tgteff"
    elif [[ "$flag" == "CU" ]] ; then
        features="VCU\n$dmax\n$tgteff"
    elif [[ "$flag" == "SU" ]] ; then
        features="VCSU\n$dmax\n$nsub\n$tgteff"
    fi
    runcard="$nenter\n$itmx\n$seed\n$npieces\nout/gen-$flag-$order-$xi-$seed.mcmule\n$features\n$pieces"
    echo -e "$runcard" | time /path/to/mcmule --gen ./user.so | tee log/gen-$flag-$order-$xi-$seed.txt
}
