#!/bin/bash
build_piece () {
    for i in $@ ; do
        echo /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/07_07_200MeV_Q2big_xi01_gen_1442/input/$i
    done
}
generate() {
    order="$1"
    flag="$2"
    xi="$3"
    seed="$4"
    tgteff=1e-3
    dmax=1.

    if [[ "$order" == 0 ]] ; then
        npieces=1
        pieces=$(build_piece mp2mp0_mu-p_S0000051157X1.00000D1.00000_ITMX020x008.0M.mcmule) #put here correct File for order 0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nenter=5000
        itmx=50
        tgteff=-1
    elif [[ "$order" == 1 ]] ; then
        npieces=2
        nenter=5000
        itmx=50
        nsub=30

        if [[ "$xi" == "01" ]] ; then
            pieces=$(build_piece \
                mp2mpNLO0_mu-p_S0000019559X0.10000D0.10000_ITMX020x008.0M.mcmule \
                mp2mpR_mu-p_S0000069163X0.10000D0.10000_ITMX020x008.0M.mcmule)     #put here correct File for order 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elif [[ "$xi" == "10" ]] ; then
            pieces=$(build_piece \
                mp2mpNLO0_blabla.mcmule \
                mp2mpR_blabla.mcmule) #put here correct File for order 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    runcard="$nenter\n$itmx\n$seed\n$npieces\nout/gen-$flag-$order-$xi-$seed.lhef\n$features\n$pieces"
    echo "===== RUNCARD ====="
    echo -e "$runcard"
    echo "==================="
    echo -e "$runcard" | /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule --gen /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/07_07_200MeV_Q2big_xi01_gen_1442/user.so
    #echo -e "$runcard" | time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule --gen /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/07_07_200MeV_Q2big_xi01_gen_1442/user.so | tee log/gen-$flag-$order-$xi-$seed.txt
}
