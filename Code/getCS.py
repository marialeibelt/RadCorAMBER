"""Extract physically normalised cross sections from the McMule event-generator branch.

The Fortran integration already uses physical alpha. The mp2mpNLO0 piece
contains Born + virtual + integrated FKS terms; mp2mpR is the FKS-subtracted
real contribution. Only their combination is an infrared-safe NLO prediction.

Do NOT interpret Delta NLO or mp2mpR as a positive event fraction.
Do NOT feed the new output into the existing event-scaling script without
reviewing the script's semantics, especially for unweighted/signed LHE events.
"""

from pymule import setup, sigma, mergefks, conv
import numpy as np
import argparse
import os
import sys


class Tee:
    """Write output to the terminal and to a file."""

    def __init__(self, filename):
        self.terminal = sys.stdout
        self.file = open(filename, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def mb_value_and_error(result):
    """McMule record.value = [central value, statistical 1-sigma error] in microbarn."""
    value, error = result.value[:2]
    return float(value) / 1000.0, abs(float(error)) / 1000.0


def get_cross_sections(homedir, outfolder, savename):
    outdir = os.path.join(homedir, outfolder)
    mcmule_outdir = os.path.join(outdir, "out")

    if not os.path.isdir(mcmule_outdir):
        raise FileNotFoundError(f"McMule output directory does not exist: {mcmule_outdir}")

    outfolder_name = os.path.basename(os.path.normpath(outfolder))
    savename_base = f"{savename}_{outfolder_name}"
    log_file = os.path.join(outdir, f"{savename_base}_output.txt")

    old_stdout = sys.stdout
    tee = Tee(log_file)
    sys.stdout = tee

    try:
        print("=" * 74)
        print("McMule cross sections:", savename_base)
        print("Input:", mcmule_outdir)
        print("=" * 74)
        setup(folder=mcmule_outdir)

        # Fortran's run_mcmule/generator sets physical alpha before integration.
        # Do NOT multiply these results by additional powers of alpha.
        lo = mergefks(sigma("mp2mp0")) * conv
        nlo0 = mergefks(sigma("mp2mpNLO0")) * conv
        real_fks = mergefks(sigma("mp2mpR")) * conv
        full = mergefks(sigma("mp2mpNLO0"), sigma("mp2mpR")) * conv

        s_lo, e_lo = mb_value_and_error(lo)
        s_nlo0, e_nlo0 = mb_value_and_error(nlo0)
        s_r, e_r = mb_value_and_error(real_fks)
        s_full, _ = mb_value_and_error(full)

        # The three pieces were integrated with independent seeds.
        # Compute uncertainty explicitly; NEVER subtract the statistical errors.
        # This formula is for independent integrations (no shared random numbers).
        e_full = np.hypot(e_nlo0, e_r)
        s_delta = s_full - s_lo
        e_delta = np.hypot(e_full, e_lo)

        delta_percent = 100.0 * s_delta / s_lo
        # delta/LO = full/LO - 1, with shared LO in numerator and denominator.
        delta_percent_error = 100.0 * np.hypot(
            e_full / s_lo, s_full * e_lo / s_lo**2
        )
        k_factor = s_full / s_lo
        k_error = delta_percent_error / 100.0

        print("All values in mb; errors are statistical (1 sigma).")
        for label, val, err in [
            ("LO", s_lo, e_lo),
            ("NLO0 (Born + virtual + integrated FKS)", s_nlo0, e_nlo0),
            ("R (FKS-subtracted; xi-dependent)", s_r, e_r),
            ("Full NLO = NLO0 + R", s_full, e_full),
            ("Delta NLO = Full - LO", s_delta, e_delta),
        ]:
            print(f"{label:43s}: {val:.12g} +/- {err:.3g} mb")

        print(f"Relative correction Delta/LO: {delta_percent:.9g} +/- "
              f"{delta_percent_error:.3g} %")
        print(f"Standard total K-factor Full/LO: {k_factor:.12g} +/- {k_error:.3g}")
        print("NOTE: R is FKS-subtracted and generally depends on xi; it is NOT")
        print("      a standalone physical hard-photon cross section.")
        print("NOTE: Delta NLO is signed. It is NOT an event mixing fraction.")

        # Diagnostic: should agree for a single xi value; if many xi values,
        # a fit/merge can produce small deviations.
        print("Check Full - (NLO0 + R) [mb]:", s_full - (s_nlo0 + s_r))

        # Keep a four-value array for convenient comparison, but use a NEW name.
        # We deliberately do NOT overwrite the old cross_sections.txt consumed
        # by any existing event-scaling script (whose semantics need rechecking).
        corrected_file = os.path.join(
            outdir, f"{savename_base}_cross_sections_corrected.txt"
        )
        cross_sections = np.array([s_lo, s_r, s_delta, s_full])
        np.savetxt(
            corrected_file,
            cross_sections,
            fmt="%.16e",
            header=(
                "Corrected McMule event-generator cross sections [mb].\n"
                "Row 1: LO (mp2mp0)\n"
                "Row 2: R (mp2mpR, FKS-subtracted; xi-dependent; diagnostic only)\n"
                "Row 3: Delta NLO = Full - LO (signed correction, NOT event fraction)\n"
                "Row 4: Full NLO = mp2mpNLO0 + mp2mpR (includes Born)\n"
                "Do not feed this file to legacy scale scripts before reviewing them."
            ),
        )
        print("Corrected cross-section array saved to:", corrected_file)

        # Save separate, labeled values and errors for slide preparation.
        error_file = os.path.join(
            outdir, f"{savename_base}_cross_sections_with_errors_corrected.txt"
        )
        np.savetxt(
            error_file,
            np.array([
                [s_lo, e_lo], [s_nlo0, e_nlo0], [s_r, e_r],
                [s_full, e_full], [s_delta, e_delta]
            ]),
            fmt="%.16e",
            header=(
                "Columns: cross_section_mb statistical_1sigma_error_mb\n"
                "Rows: LO; NLO0 (Born inclusive); R (FKS/xi-dependent); "
                "Full NLO; Delta NLO"
            )
        )
        print("Cross sections and errors saved to:", error_file)
    finally:
        sys.stdout = old_stdout
        tee.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate corrected McMule cross sections")
    parser.add_argument("--homedir", required=True)
    parser.add_argument("--outfolder", required=True)
    parser.add_argument("--savename", required=True)
    args = parser.parse_args()
    get_cross_sections(args.homedir, args.outfolder, args.savename)
