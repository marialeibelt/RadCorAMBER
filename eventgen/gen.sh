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
    dmax=0.0045 #put here correct value for dmax!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if [[ "$order" == 0 ]] ; then
        npieces=1
        pieces=$(build_piece mp2mp0_blabla.vegas) #put here correct File for order 0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                mp2mpNLO0_blabla.vegas \
                mp2mpR_blabla.vegas)     #put here correct File for order 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elif [[ "$xi" == "10" ]] ; then
            pieces=$(build_piece \
                mp2mpNLO0_blabla.vegas \
                mp2mpR_blabla.vegas) #put here correct File for order 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    runcard="$nenter\n$itmx\n$seed\n$npieces\nout/gen-$flag-$order-$xi-$seed.mme\n$features\n$pieces"
    echo -e "$runcard" | time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-release/mcmule --gen ./user.so | tee log/gen-$flag-$order-$xi-$seed.txt
}
