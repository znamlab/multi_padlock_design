# padlock design pipeline for multiplexed assay with multiple probes
# per target in cDNA-based expression profiling
# Xiaoyan, 2017
import sys
import traceback

import multi_padlock_design.config as config
from multi_padlock_design.blast import parblast, readblast
from multi_padlock_design.io import checkinput, createoutput, formatrefseq
from multi_padlock_design.select_probes import (
    distributeprobes,
    finalizeprobes,
    screenseq,
)

MODE = "STARBARseq"  #'STARBARseq' or 'BARseq'


def main() -> int:
    """Run the interactive probe design pipeline."""

    try:
        # get keyboard inputs and prepare sequences

        # design_pars : (species, int(armlen), int(interval), int(t1), int(t2),
        # number of probes per gene, total length)
        # outpars : (outdir, outdir_temp, input_type)
        # genepars : (genes, linkers, headers, variants)
        # designinput : (basepos, headers_wpos, sequences, variants_matching_sequence)

        # get keyboard inputs and prepare sequences
        designpars, outpars, genepars, designinput = checkinput.getdesigninput()
        # fmt = checkinput.checkformat(genepars[2])

        # Tm screening
        Tm, siteChopped = screenseq.thresholdtm(
            designinput[1],
            designinput[2],
            outpars[1],
            designpars,  # armlen used, but not important, only total length
        )

        # Make blast database if it doesn't exist
        formatrefseq.blastdb(designpars[0])
        # Run blast in parallel
        parblast.continueblast(siteChopped, designinput[1], outpars[1], designpars)

        # Find specific targets
        siteCandidates, notMapped = readblast.getcandidates(
            siteChopped,
            designinput[1],
            outpars,
            designpars[1],
            designinput[3],
            config.specificity_by_tm,
            designpars[6],  # total length
        )  # used in readblastout
        print(f"Not mapped genes after candidate assignment: {notMapped}")
        createoutput.writetargetfile(
            designinput,
            siteCandidates,
            Tm,
            designpars[1],
            outpars,
            "3.AllSpecificTargets_",
            designpars[6],  # total length
        )

        # non-overlapping candidates
        targets, targetpos, mapTmlist = distributeprobes.asmanyprobes(
            siteCandidates, siteChopped, designinput[2], designpars
        )

        # correct positions
        targets, targetpos, notMapped, Tm = finalizeprobes.correctpos(
            designinput[0], targets, targetpos, notMapped, mapTmlist, Tm, siteChopped
        )

        # write genes with no candidates
        createoutput.emptyentries(targets, genepars[2], outpars)

        # fill up linker sequence and write
        probes = finalizeprobes.assembleprobes(targets, genepars, designpars[1])
        createoutput.writeprobefile(
            genepars[0],
            genepars[2],
            probes,
            Tm,
            targetpos,
            targets,
            outpars,
            designpars[1],
            "4.NonOverlappingProbes_",
        )

        # remove targets cannot be found in database
        finallist = finalizeprobes.removeunmapped(
            notMapped, targetpos, genepars[2], targets, Tm, probes
        )
        createoutput.writeprobefile(
            genepars[0],
            genepars[2],
            finallist[0],
            finallist[1],
            finallist[2],
            finallist[3],
            outpars,
            designpars[1],
            "5.ProbesDBMappable_",
        )

        # prioritize sequences without homopolymers and randomly select the fixed number
        #  of probes per gene (if applicable)
        if len(designpars[5]):
            sublist = finalizeprobes.selectprobes(
                int(designpars[5]), finallist, genepars[2], designpars[1], outpars
            )
            createoutput.writeprobefile(
                genepars[0],
                genepars[2],
                sublist[0],
                sublist[1],
                sublist[2],
                sublist[3],
                outpars,
                designpars[1],
                "6.ProbesRandomSubsetN=" + designpars[5] + "_",
                regions=sublist[4],
            )

        print("All finished!")
        return 0

    except Exception:
        print(sys.exc_info()[0])
        print(traceback.format_exc())
        return 1

    finally:
        if sys.stdin is not None and sys.stdin.isatty():
            print("Press Enter to continue ...")
            input()


if __name__ == "__main__":
    sys.exit(main())
