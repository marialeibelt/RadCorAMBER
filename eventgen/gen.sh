#!/bin/bash
basefolder="/nfs/momos/user/mleibelt"   # "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"

# ---------------change this each time you want to run a new analysis---------------
anafolder="05_08_evtgen_25_08" # "08_07_200MeV_Q2big_xi01_gentest_1444"

mcmule0="mp2mp0_mu-p_S0000077682X1.00000D1.00000_ITMX020x080.0M.mcmule" # "mp2mp0_mu-p_S0000085728X1.00000D1.00000_ITMX020x008.0M.mcmule"

mcmuleNLO0_xi01="mp2mpNLO0_mu-p_S0000075957X0.10000D0.10000_ITMX020x080.0M.mcmule" # "mp2mpNLO0_mu-p_S0000044063X0.10000D0.10000_ITMX020x008.0M.mcmule"
mcmuleR_xi01="mp2mpR_mu-p_S0000023375X0.10000D0.10000_ITMX020x080.0M.mcmule" # "mp2mpR_mu-p_S0000040361X0.10000D0.10000_ITMX020x008.0M.mcmule"

mcmuleNLO0_xi10="mp2mpNLO0_blabla.mcmule"
mcmuleR_xi10="mp2mpR_blabla.mcmule"
# ----------------------------------------------------------------------------------

build_piece () {
    for i in $@ ; do
        echo "$basefolder/$anafolder/input/$i"
    done
}
mkdir -p "$basefolder/$anafolder/log"
generate() {
    order="$1"
    flag="$2"
    xi="$3"
    seed="$4"
    tgteff=1e-1
    dmax=1.

    if [[ "$order" == 0 ]] ; then
        npieces=1
        pieces=$(build_piece "$mcmule0")
        nenter=500
        itmx=12
        tgteff=-1

    elif [[ "$order" == 1 ]] ; then
        npieces=2
        nenter=500
        itmx=12
        nsub=30

        if [[ "$xi" == "01" ]] ; then
            pieces=$(build_piece \
                "$mcmuleNLO0_xi01" \
                "$mcmuleR_xi01")

        elif [[ "$xi" == "10" ]] ; then
            pieces=$(build_piece \
                "$mcmuleNLO0_xi10" \
                "$mcmuleR_xi10")
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
    runcard="$nenter\n$itmx\n$seed\n$npieces\nout/gen-$flag-$order-$xi-$seed.lhe\n$features\n$pieces"

    echo "===== RUNNING <3 ====="
    echo -e "$runcard" | time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule --gen "$basefolder/$anafolder/user.so" | tee "$basefolder/$anafolder/log/gen-$flag-$order-$xi-$seed.txt"
}