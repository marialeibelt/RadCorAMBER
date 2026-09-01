from pymule import *
import numpy as np
import argparse
import os
import sys


class Tee:
    """Writes output to both terminal and a file."""

    def __init__(self, filename):
        self.terminal = sys.stdout
        self.file = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def get_cross_sections(homedir, outfolder, savename):

    # ==========================================
    # Paths
    # ==========================================

    outdir = os.path.join(homedir, outfolder)
    mcmule_outdir = os.path.join(outdir, "out")

    if not os.path.isdir(mcmule_outdir):
        raise FileNotFoundError(
            f"McMule output directory does not exist:\n{mcmule_outdir}"
        )

    # Only use folder name for output file names
    outfolder_name = os.path.basename(os.path.normpath(outfolder))

    savename_base = f"{savename}_{outfolder_name}"

    # ==========================================
    # Redirect stdout
    # ==========================================

    log_file = os.path.join(
        outdir,
        f"{savename_base}_output.txt"
    )

    old_stdout = sys.stdout
    tee = Tee(log_file)
    sys.stdout = tee

    try:

        print("=======================================================================")
        print("Analysed file:", savename_base)
        print("=======================================================================")
        print("McMule output directory:", mcmule_outdir)
        print()

        # ==========================================
        # Read McMule results
        # ==========================================

        setup(folder=mcmule_outdir)

        lo = (
            mergefks(
                sigma("mp2mp0")
            )
            * alpha**2
            * conv
        )

        nlo = (
            mergefks(
                sigma("mp2mpR"),
                sigma("mp2mpNLO0")
            )
            * alpha**3
            * conv
        )

        onlyR = (
            mergefks(
                sigma("mp2mpR")
            )
            * alpha**3
            * conv
        )

        full = lo + nlo

        # ==========================================
        # Convert to mb
        # ==========================================

        cs_lo = lo.value[0] / 1000
        cs_onlyR = onlyR.value[0] / 1000
        cs_nlo = nlo.value[0] / 1000
        cs_full = full.value[0] / 1000

        # ==========================================
        # Output
        # ==========================================

        print("LO cross section:     ", cs_lo, "mb")
        print("NLO cross section:    ", cs_nlo, "mb")
        print("Full cross section:   ", cs_full, "mb")
        print("Only R cross section: ", cs_onlyR, "mb")

        cs_array = np.array([
            cs_lo,
            cs_onlyR,
            cs_nlo,
            cs_full
        ])

        print()
        print(
            "Cross section array "
            "(LO, onlyR, NLO, full):",
            cs_array
        )

        # ==========================================
        # Save array
        # ==========================================

        array_file = os.path.join(
            outdir,
            f"{savename_base}_cross_sections.txt"
        )

        np.savetxt(
            array_file,
            cs_array,
            header="1. LO  2. onlyR  3. NLO  4. full"
        )

        print()
        print("Cross section array saved to:")
        print(array_file)

    finally:
        sys.stdout = old_stdout
        tee.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate McMule cross sections"
    )

    parser.add_argument(
        "--homedir",
        required=True
    )

    parser.add_argument(
        "--outfolder",
        required=True
    )

    parser.add_argument(
        "--savename",
        required=True
    )

    args = parser.parse_args()

    get_cross_sections(
        homedir=args.homedir,
        outfolder=args.outfolder,
        savename=args.savename
    )