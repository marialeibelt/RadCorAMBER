#!/bin/bash

# Path of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Read configuration
source "$SCRIPT_DIR/input_evtgen_mcmule.conf"


build_piece () {

    for i in "$@"; do
        echo "$basefolder/$anafolder/input/$i"
    done
}


mkdir -p "$basefolder/$anafolder/log"
mkdir -p "$basefolder/$anafolder/out"


generate() {

    order="$1"
    flag="$2"
    xi="$3"
    seed="$4"


    # ============================================================
    # LO
    # ============================================================

    if [[ "$order" == "0" ]]; then

        npieces=1

        nenter="$nenter_LO"
        itmx="$itmx_LO"
        tgteff="$tgteff_LO"

        if [[ "$xi" == "01" ]]; then

            pieces=$(build_piece \
                "$mcmule0_xi01"
            )

        elif [[ "$xi" == "10" ]]; then

            pieces=$(build_piece \
                "$mcmule0_xi10"
            )

        else

            echo "ERROR: Unknown xi = $xi"
            return 1

        fi


    # ============================================================
    # NLO
    # ============================================================

    elif [[ "$order" == "1" ]]; then

        npieces=2

        nenter="$nenter_NLO"
        itmx="$itmx_NLO"
        tgteff="$tgteff_NLO"
        nsub="$nsub_NLO"

        if [[ "$xi" == "01" ]]; then

            pieces=$(build_piece \
                "$mcmuleNLO0_xi01" \
                "$mcmuleR_xi01"
            )

        elif [[ "$xi" == "10" ]]; then

            pieces=$(build_piece \
                "$mcmuleNLO0_xi10" \
                "$mcmuleR_xi10"
            )

        else

            echo "ERROR: Unknown xi = $xi"
            return 1

        fi

    else

        echo "ERROR: Unknown order = $order"
        return 1

    fi


    # ============================================================
    # Generator flags
    # ============================================================

    if [[ "$flag" == "N" ]]; then

        features="V"

    elif [[ "$flag" == "C" ]]; then

        features="VC\n$dmax"

    elif [[ "$flag" == "S" ]]; then

        features="VCS\n$dmax\n$nsub"

    elif [[ "$flag" == "U" ]]; then

        features="VU\n$tgteff"

    elif [[ "$flag" == "CU" ]]; then

        features="VCU\n$dmax\n$tgteff"

    elif [[ "$flag" == "SU" ]]; then

        features="VCSU\n$dmax\n$nsub\n$tgteff"

    else

        echo "ERROR: Unknown flag = $flag"
        return 1

    fi


    # ============================================================
    # Runcard
    # ============================================================

    runcard="$nenter
$itmx
$seed
$npieces
out/gen-$flag-$order-$xi-$seed.lhe
$features
$pieces"


    # ============================================================
    # Run generator
    # ============================================================

    echo -e "$runcard" | \
        time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule \
        --gen "$basefolder/$anafolder/user.so" \
        | tee "$basefolder/$anafolder/log/gen-$flag-$order-$xi-$seed.txt"
}