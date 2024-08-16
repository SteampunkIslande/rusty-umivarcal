import math
import os
import re
import statistics
import sys
import time
from collections import OrderedDict
from multiprocessing import Pool

import msgpack

# import local modules
from functions import *
from pyfasta import Fasta

# from ParseBED import ParseBED
# from FilterPositions import FilterPositions
# from CallVariants import CallVariants
# from Phase import Phase
# from func import Exit_dev
# from func import ConcatenatePileups
# from func import SortPileup
# from func import Add_Depth_Noise_Ref_HP
# from func import PrintTime
# from func import GetTotalLines
# from func import GetTotalVariants
# from func import *

# from TreatReads import TreatReads
# from Output import Output
# from CalculateStats import CalculateStats
# from MergeSubPileups import MergeSubPileups
# from ComputeCoverageStats import ComputeCoverageStats


def Call(config, FUNC_PATH):

    # load and define variables from the config
    INPUT = config["input"]
    BED = config["bed"]
    FASTA = config["fasta"]
    MIN_BASE_QUALITY = int(config["min_base_quality"])
    MIN_MAPPING_QUALITY = int(config["min_mapping_quality"])
    MIN_READ_QUALITY = int(config["min_read_quality"])
    MIN_VARIANT_UMI = int(config["min_variant_umi"])
    MIN_VARIANT_DEPTH = int(config["min_variant_depth"])
    STRAND_BIAS_METHOD = str(config["strand_bias_method"])
    MAX_STRAND_BIAS = float(config["max_strand_bias"])
    PILEUP = config["pileup"]
    REBUILD = False if os.path.isfile(PILEUP) else True
    OUTPUT = config["output"]
    CORES = config["cores"]
    DEFAULT_CORES = config["default_cores"]
    ALPHA = float(config["alpha"])
    MAX_HP_LENGTH = int(config["max_hp_length"])
    gVCF = config["gvcf"]
    KEEP_PILEUP = config["keep_pileup"]
    BLACK_LIST = config["black_list"]
    WHITE_LIST = config["white_list"]
    MIN_PHASE_UMI = int(config["min_phase_umi"])
    MIN_PHASE_VAF_RATIO = float(config["min_phase_vaf_ratio"])
    MAX_PHASE_DISTANCE = int(config["max_phase_distance"])
    COMPUTE_COVERAGE_STATS = config["compute_coverage_stats"]
    VERSION = config["version"]
    LAST_UPDATE = config["lastUpdate"]

    # print parameters in the console
    PrintTime("green", "\t\tINPUT file     : " + INPUT)
    PrintTime("green", "\t\tBED file       : " + BED)
    PrintTime("green", "\t\tFASTA file     : " + FASTA)

    if PILEUP != "None":
        PrintTime("green", "\t\tPILEUP file    : " + PILEUP)
    PrintTime("green", "\t\tOutput         : " + OUTPUT)

    PrintTime("green", "\t\tmin_base_quality       : " + str(MIN_BASE_QUALITY))
    PrintTime("green", "\t\tmin_read_quality       : " + str(MIN_READ_QUALITY))
    PrintTime("green", "\t\tmin_mapping_quality    : " + str(MIN_MAPPING_QUALITY))
    PrintTime("green", "\t\tmin_variant_umi        : " + str(MIN_VARIANT_UMI))
    PrintTime("green", "\t\tmin_variant_depth      : " + str(MIN_VARIANT_DEPTH))
    PrintTime("green", "\t\tstrand_bias_method     : " + str(STRAND_BIAS_METHOD))
    PrintTime("green", "\t\tmax_strand_bias        : " + str(MAX_STRAND_BIAS))
    PrintTime("green", "\t\tmax_hp_length          : " + str(MAX_HP_LENGTH))
    PrintTime("green", "\t\talpha                  : " + str(ALPHA))
    PrintTime("green", "\t\tmin_phase_umi          : " + str(MIN_PHASE_UMI))
    PrintTime("green", "\t\tmin_phase_vaf_ratio    : " + str(MIN_PHASE_VAF_RATIO))
    PrintTime("green", "\t\tmax_phase_distance     : " + str(MAX_PHASE_DISTANCE))
    if gVCF:
        PrintTime(
            "green", "\t\tgVCF                   : " + str(gVCF) + " (Experimental)"
        )
    else:
        PrintTime("green", "\t\tgVCF                   : " + str(gVCF))

    if DEFAULT_CORES:
        PrintTime("green", "\t\tcores                  : " + str(CORES) + " (default)")
    else:
        PrintTime("green", "\t\tcores                  : " + str(CORES))

    PrintTime("green", "\t\tkeep_pileup            : " + str(KEEP_PILEUP))
    PrintTime("green", "\t\tcompute_coverage_stats : " + str(COMPUTE_COVERAGE_STATS))
    PrintTime("green", "\t\tblack_list             : " + str(BLACK_LIST))
    PrintTime("green", "\t\twhite_list             : " + str(WHITE_LIST) + "\n")

    PrintTime("green", "\t\tVERSION                : " + VERSION)

    PrintTime("console", "\tDone\n")

    # make dir if outdir doesn't exist
    try:
        os.mkdir(OUTPUT)
    except:
        pass

    # load the reference genome file
    f = Fasta(FASTA)

    # if input is bam => launch samtools view command
    # to convert it to sam
    if ".bam" in INPUT and ".sam" not in INPUT:

        print("\n")
        PrintTime("console", "\tConverting BAM to SAM...")

        SAM = BAMtoSAM(INPUT)

        PrintTime("console", "\tDone")

    else:
        # else => sam = input
        SAM = INPUT

    nReads = GetTotalLines(SAM)
    validReads = 0

    # if a pileup is not given, the pileup has to be build
    if REBUILD:

        print("\n")
        PrintTime("console", "\tAnalyzing BED...")

        # build the empty pileup
        pileup = ParseBED(BED)

        pileupLen = GetPileupLength(pileup)
        # print(pileupLen); exit()

        if CORES > 1:

            # if pileup length < 1M positions => using multiple cores is advantageous
            # if 1M < pileup length <= 2M => a minimum of 5M reads is required to see performance gains
            # if 2M < pileup length <= 5M => a minimum of 10M reads is required to see performance gains
            # if 5M < pileup length <= 10M => a minimum of 50M reads is required to see any performance gains
            # if pileup length > 10M, theoretically, 100M reads should be analyzed faster but the memory usage will skyrocket so automatically switch to 1 core only.
            if (
                (pileupLen > 1000000 and pileupLen <= 2000000 and nReads < 5000000)
                or (pileupLen > 2000000 and pileupLen <= 5000000 and nReads < 20000000)
                or (pileupLen > 5000000 and pileupLen <= 10000000 and nReads < 50000000)
                or (pileupLen > 10000000)
            ):

                PrintTime(
                    "warning",
                    "\t\tWarning: Using more cores to analyze the provided data will not result in any significant performance gains!\n\t\t\t\tLaunching UMI-VarCal on one core only...\n",
                )
                CORES = 1

        # dump the pileup object
        with open(OUTPUT + "/.pileup_ini", "wb") as handle:
            msgpack.pack(pileup, handle, encoding="utf-8")

        PrintTime("console", "\tDone\n")

        if KEEP_PILEUP:
            PrintTime(
                "warning",
                "\tWarning: keep_pileup parameter is set to True. This can affect your memory consumption.\n\t\t\tIf your system runs out of memory, please try to run the analysis with --keep_pileup False instead.\n",
            )
        PrintTime("console", "\tBuilding Pileup...")

        # if multiple cores are used
        # wait until all processes are finished == until all files are generated
        if CORES > 1:

            #############################################################################################
            ###########################                                       ###########################
            ########################### PARALLELIZED CODE START : TREAT READS ###########################
            ###########################                                       ###########################
            #############################################################################################

            # split the SAM files into equal sub files
            nReads_split = int(nReads / float(CORES))

            # preprocess reads
            # if more then one core is to be used, separate the input into subfiles
            subFiles = PreprocessReads(SAM, CORES)

            # parallelization block
            p = Pool(processes=int(CORES))
            tmp = []
            for subFile in subFiles:
                # get the values from each core and add return values to the tmp array
                tmp.append(
                    p.apply_async(
                        TreatReads,
                        (
                            subFile,
                            GetTotalLines(subFile),
                            OUTPUT,
                            MIN_BASE_QUALITY,
                            MIN_READ_QUALITY,
                            MIN_MAPPING_QUALITY,
                        ),
                    )
                )
            p.close()
            p.join()

            # parse the tmp array and get pileups and valid reads from each core
            pileups = []
            subValidReads = []
            for x in tmp:
                pileups.append(x.get()[0])
                subValidReads.append(x.get()[1])

            # total valid reads = the sum of valid reads from each core
            validReads = sum(subValidReads)

            # merge sub pileups to obtain whole pileup
            pileup = MergeSubPileups(pileup, pileups, subFiles, OUTPUT)

            #########################################################################################
            ###########################                                     #########################
            ########################### PARALLELIZED CODE END : TREAT READS #########################
            ###########################                                     #########################
            #########################################################################################

        else:

            # if only one core is to used, launch the function from here since no need to merge
            value = TreatReads(
                SAM,
                nReads,
                OUTPUT,
                MIN_BASE_QUALITY,
                MIN_READ_QUALITY,
                MIN_MAPPING_QUALITY,
            )
            pileup = value[0]
            validReads = value[1]

        print("\n")
        PrintTime("console", "\tDone")

        print("\n")
        PrintTime("console", "\tEstimating Noise in Reads...")

        # add depth to pileup
        # add variant error noise at each position
        # add reference bases in the dictionnary
        # add homopolymers infos

        pileup = Add_Depth_Noise_Ref_HP(pileup, FASTA)

        # rebuild to SAM original name
        SAM = SAM.replace("_reordered.sam", ".sam")

        # dump pileup in msgpack object
        if KEEP_PILEUP:
            with open(
                OUTPUT + "/" + SAM.replace(".sam", ".pileup").split("/")[-1], "wb"
            ) as handle:
                msgpack.pack(pileup, handle, encoding="utf-8")

        # remove tmp empty pileup file
        try:
            os.remove(OUTPUT + "/.pileup_ini")
        except:
            pass

        print("\n")
        PrintTime("console", "\tDone")

    else:

        print("\n")
        PrintTime("console", "\tLoading Pileup...")

        # load pileup from msgpack object
        with open(PILEUP, "rb") as handle:
            pileup = msgpack.unpack(
                handle, encoding="utf-8", max_buffer_size=100 * 1024 * 1024 * 1024
            )
            pileup = SortPileup(pileup)

        PrintTime("console", "\tDone")

    ### Compute UMI coverage statistics per region & position
    if COMPUTE_COVERAGE_STATS:
        ComputeCoverageStats(pileup, BED, OUTPUT, SAM)

    pileupLen = GetPileupLength(pileup)

    # print(pileup)
    full_pileup = CopyPileup(pileup)

    ### parse black and white lists if provided
    if BLACK_LIST != "None" or WHITE_LIST != "None":
        print("\n")
        PrintTime("console", "\tParsing list(s)...")

        value = ParseLists(BLACK_LIST, WHITE_LIST)
        BLACK_LIST = value[0]
        WHITE_LIST = value[1]

        PrintTime("console", "\tDone")
    else:
        BLACK_LIST = []
        WHITE_LIST = {"pos": [], "var": []}

    ### Poisson modeling to filter positions
    result = FilterPositions(pileup, ALPHA, WHITE_LIST)
    pileup = result[0]
    potential = result[1]
    foundCandidates = result[2]

    if foundCandidates:

        ### call final variants
        value = CallVariants(
            pileup,
            f,
            STRAND_BIAS_METHOD,
            MAX_STRAND_BIAS,
            MIN_VARIANT_UMI,
            MIN_VARIANT_DEPTH,
            MAX_HP_LENGTH,
            BLACK_LIST,
            WHITE_LIST,
        )
        finalVariants = value[0]
        phasedVariants = value[1]

        ### call phased variants
        phasedVariants = Phase(
            phasedVariants, MIN_PHASE_UMI, MIN_PHASE_VAF_RATIO, MAX_PHASE_DISTANCE
        )

        ### Writing results to VCF
        final = Output(
            full_pileup,
            pileup,
            finalVariants,
            phasedVariants,
            INPUT,
            SAM,
            BED,
            FASTA,
            OUTPUT,
            MIN_BASE_QUALITY,
            MIN_READ_QUALITY,
            MIN_MAPPING_QUALITY,
            MIN_VARIANT_UMI,
            STRAND_BIAS_METHOD,
            MAX_STRAND_BIAS,
            CORES,
            ALPHA,
            MAX_HP_LENGTH,
            gVCF,
            VERSION,
        )

        # calculate and display stats
        CalculateStats(pileup, potential, final, pileupLen, nReads, validReads, REBUILD)
    else:
        print("\n")
        message = "No candidate positions were found !\n"
        PrintTime("error", "\t" + message)

        ### Writing results to VCF (even if no variants found)
        finalVariants = {}
        phasedVariants = {}
        final = Output(
            full_pileup,
            pileup,
            finalVariants,
            phasedVariants,
            INPUT,
            SAM,
            BED,
            FASTA,
            OUTPUT,
            MIN_BASE_QUALITY,
            MIN_READ_QUALITY,
            MIN_MAPPING_QUALITY,
            MIN_VARIANT_UMI,
            STRAND_BIAS_METHOD,
            MAX_STRAND_BIAS,
            CORES,
            ALPHA,
            MAX_HP_LENGTH,
            gVCF,
            VERSION,
        )

        print("\r")
        PrintTime("console", "\tCalculating statistics...")

        # print out stats to console
        message = "Candidate Positions: 0"
        PrintTime("green", "\t\t" + message)
        message = "Final Variants: 0"
        PrintTime("green", "\t\t" + message)
