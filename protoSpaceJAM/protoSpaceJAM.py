import datetime
import io
import linecache
import logging
import os.path
import pickle
import json
import sys
import traceback
import time
import re
import ssl
import pandas as pd
import argparse
from types import SimpleNamespace
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.GenBank import Record


from protoSpaceJAM.util.utils import MyParser, ColoredLogger, read_pickle_files, cal_elapsed_time, get_gRNAs,get_gRNAs_target_coordinate, get_gRNAs_near_loc, get_chopchop_raw_results, convert_chopchop_raw_to_psj, get_start_stop_loc, get_end_pos_of_ATG, get_start_pos_of_stop, \
    get_HDR_template #uncomment this for pip installation

# from util.utils import MyParser, ColoredLogger, read_pickle_files, cal_elapsed_time, get_gRNAs, get_gRNAs_target_coordinate, \
#     get_HDR_template

def parse_args(test_mode=False):
    parser = MyParser(description="protoSpaceJAM: perfectionist CRISPR knock-in design at scale\n",
                      formatter_class=argparse.RawTextHelpFormatter)  # Ensures newlines are preserved
    IO = parser.add_argument_group('input/output')
    IO.add_argument(
        "--path2csv",
        default=os.path.join("input","test_input.csv"),
        type=str,
        help="Path to a csv file containing the input knock-in sites, see input/test_input.csv for an example\n *required columns*: 'Ensembl_ID'\n coordinate-based rows use: 'Chromosome','Coordinate'\n transcript-terminus rows use: 'Target_terminus'\n optional region-targeting columns: 'Preferred_region','Preferred_exon','Preferred_anchor','Preferred_offset','Gene_Name'",
        metavar="<PATH_TO_CSV>",
    )
    IO.add_argument(
        "--input_mode",
        default="csv",
        type=str,
        help="Input mode: csv or direct (default: csv). direct allows passing one target from CLI args without preparing a CSV file.",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_ensembl_id",
        default="",
        type=str,
        help="[input_mode=direct] ENST ID (optional for coordinate-based design)",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_target_terminus",
        default="ALL",
        type=str,
        help="[input_mode=direct] N, C, or ALL for ENST-based design (default: ALL)",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_chromosome",
        default="",
        type=str,
        help="[input_mode=direct] Chromosome for coordinate-based design (for example: 7 or chr7)",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_coordinate",
        default="",
        type=str,
        help="[input_mode=direct] Genomic coordinate for coordinate-based design",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_gene_name",
        default="",
        type=str,
        help="[input_mode=direct] Gene symbol/name for CHOPCHOP gene-based web query (for example: TRAC).",
        metavar="<string>",
    )
    IO.add_argument(
        "--direct_entry",
        default="1",
        type=str,
        help="[input_mode=direct] Optional entry label to write into outputs (default: 1)",
        metavar="<string>",
    )
    IO.add_argument(
        "--outdir",
        default=os.path.join("output","test"),
        type=str,
        metavar = "<PATH_TO_OUTPUT_DIRECTORY>",
        help="Path to the output directory"
    )
    IO.add_argument(
        "--guides_csv",
        default="",
        type=str,
        help="Optional path to a precomputed guide CSV in protoSpaceJAM-compatible format. "
             "If provided, guide discovery is skipped and guides are read from this file for HDR generation.",
        metavar="<PATH_TO_GUIDES_CSV>",
    )
    genome = parser.add_argument_group('genome')
    genome.add_argument(
        "--genome_ver",
        default="GRCh38",
        type=str,
        help="Genome and version to use, possible values are GRCh38, GRCm39, GRCz11, and mRatBN7.2",
        metavar="<string>",
    )
    genome.add_argument(
        "--annotation_source",
        default="auto",
        type=str,
        help="Annotation source: auto, local, or ensembl. auto uses local pickles when present; otherwise Ensembl fallback for guides_csv workflows.",
        metavar="<string>",
    )
    genome.add_argument(
        "--ensembl_cache_dir",
        default="",
        type=str,
        help="Optional directory for transcript-level Ensembl annotation cache.",
        metavar="<string>",
    )
    genome.add_argument(
        "--cache_ensembl_annotation",
        default=False,
        action="store_true",
        help="If set, write Ensembl-derived transcript annotation bundles to --ensembl_cache_dir for reuse.",
    )
    gRNA = parser.add_argument_group('gRNA')
    gRNA.add_argument(
        "--pam",
        default="NGG",
        type=str,
        help="PAM sequence (default: NGG)",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--guide_source",
        default="chopchop",
        type=str,
        help="Guide source: chopchop or precomputed (default: chopchop)",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--chopchop_cmd_template",
        default="",
        type=str,
        help="Optional local command template used to run CHOPCHOP. Must produce a tabular guide file at {output_file}. "
             "Template variables: {fasta},{output_file},{output_dir},{pam},{genome},{chrom},{pos},{window_start},{window_end}",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--chopchop_web_base_url",
        default="https://chopchop.cbu.uib.no",
        type=str,
        help="Base URL for CHOPCHOP web service (default: https://chopchop.cbu.uib.no)",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--chopchop_web_timeout",
        default=180,
        type=int,
        help="Timeout in seconds for CHOPCHOP web submit+poll workflow (default: 180)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--chopchop_web_poll_interval",
        default=2,
        type=int,
        help="Polling interval in seconds for CHOPCHOP web results (default: 2)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--chopchop_web_payload_json",
        default="",
        type=str,
        help="Optional raw JSON payload template for CHOPCHOP web submit. If omitted, a default payload is generated automatically.",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--chopchop_debug_dump",
        default=False,
        action="store_true",
        help="Write CHOPCHOP debug artifacts (input sequence window and submit payload) to <outdir>/chopchop_debug.",
    )
    gRNA.add_argument(
        "--psj_rank_chopchop",
        default=False,
        action="store_true",
        help="If set, re-rank CHOPCHOP guides with protoSpaceJAM scoring. Default keeps CHOPCHOP ranking order.",
    )
    gRNA.add_argument(
        "--chopchop_window_padding",
        default=80,
        type=int,
        help="Extra bp padding on each side of the target window fetched from Ensembl for CHOPCHOP (default: 80)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--ensembl_timeout",
        default=30,
        type=int,
        help="Timeout in seconds for Ensembl REST requests (default: 30)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--num_gRNA_per_design",
        default=1,
        type=int,
        help="Number of gRNAs to return per site (default: 1)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--max_cut2ins_dist",
        default=50,
        type=int,
        help="Maximum absolute cut-to-insert distance allowed for guide selection (default: 50)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--no_regulatory_penalty",
        default=False,
        action="store_true",
        help="Turn off penalty for gRNAs cutting in UTRs or near splice junctions, default: penalty on",
    )
    gRNA.add_argument(
        "--alpha1",
        default=1.0,
        type=float,
        help="raise the specificity weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--alpha2",
        default=1.0,
        type=float,
        help="raise the insert dist. weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--alpha3",
        default=1.0,
        type=float,
        help="raise the position weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--specificity_backend",
        default="default",
        type=str,
        help="Specificity scoring backend for CHOPCHOP-derived guides: default, chopchop_proxy, or crispor (default: default)",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--chopchop_proxy_k",
        default=600.0,
        type=float,
        help="Half-penalty point for the CHOPCHOP proxy weight: weight = 1/(1 + (penalty/K)^n) (default: 600)",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--chopchop_proxy_n",
        default=1.0,
        type=float,
        help="Steepness exponent n for the CHOPCHOP proxy weight (default: 1.0)",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--crispor_cmd_template",
        default="",
        type=str,
        help="Optional local command template used to score guides with CRISPOR. Must write a tabular file at {output_file}. Template variables: {input_file},{output_file},{output_dir},{genome},{pam}",
        metavar="<string>",
    )
    payload = parser.add_argument_group('payload')
    payload.add_argument(
        "--payload_type",
        default="insertion",
        type=str,
        help="Define the payload type, possible values are 'insertion', 'SNP'\n"
        "if payload_type is 'insertion', the payload sequence is the insertion sequence. payload sequence should be defined with --payload, or --Npayload, --Cpayload, --POSpayload, --Tag, --Linker.\n"
        "if payload_type is 'SNP', the payload sequence is the SNP sequence. The payload sequence will replace the sequence starting at the specified coordinate defined by 'Chromosome' and 'Coordinate' columns of the input csv file. The payload sequence should be defined with --SNPpayload. ",
        metavar="<string>",
    )
    payload.add_argument(
        "--payload",
        default="",
        type=str,
        help="Define the payload sequence for every site, regardless of terminus or coordinates, overrides --Npayload, --Cpayload, POSpayload, --Tag, --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--payload_file",
        default="",
        type=str,
        help="Optional path to a text/FASTA file containing payload sequence for all designs (same behavior as --payload).",
        metavar="<PATH_TO_FILE>",
    )
    payload.add_argument(
        "--Npayload",
        default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        type=str,
        help="Payload sequence to use at the N terminus (default: mNG11 + XTEN80): ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--Npayload_file",
        default="",
        type=str,
        help="Optional path to a text/FASTA file containing N-terminus payload sequence (same behavior as --Npayload).",
        metavar="<PATH_TO_FILE>",
    )
    payload.add_argument(
        "--Cpayload",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="Payload sequence to use at the C terminus (default: XTEN80 + mNG11): GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--Cpayload_file",
        default="",
        type=str,
        help="Optional path to a text/FASTA file containing C-terminus payload sequence (same behavior as --Cpayload).",
        metavar="<PATH_TO_FILE>",
    )
    payload.add_argument(
        "--POSpayload",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="Payload sequence to use at the specific genomic coordinates (default: XTEN80 + mNG11): GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--POSpayload_file",
        default="",
        type=str,
        help="Optional path to a text/FASTA file containing coordinate-mode payload sequence (same behavior as --POSpayload).",
        metavar="<PATH_TO_FILE>",
    )
    payload.add_argument(
        "--Tag",
        default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="default is the mNG11 tag",
        metavar="<string>",
    )
    payload.add_argument(
        "--Linker",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        type=str,
        help="default is the XTEN80 linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--SNPpayload",
        default="",
        type=str,
        help="Payload sequence to use for SNPs",
    )
    payload.add_argument(
        "--SNPpayload_file",
        default="",
        type=str,
        help="Optional path to a text/FASTA file containing SNP payload sequence (same behavior as --SNPpayload).",
        metavar="<PATH_TO_FILE>",
    )

    donor = parser.add_argument_group('donor')
    donor.add_argument(
        "--Donor_type",
        default="ssODN",
        help="Set the type of donor, possible values are ssODN and dsDNA (default: ssODN)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--HA_len",
        default=500,
        help="[dsDNA] Length of the desired homology arm on each side (default: 500)",
        type=int,
        metavar="<integer>",
    )
    donor.add_argument(
        "--Strand_choice",
        default="auto",
        help="[ssODN] Strand choice of ssoODN, Possible values are 'auto', 'TargetStrand', 'NonTargetStrand', 'CodingStrand' and 'NonCodingStrand'",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--ssODN_max_size",
        type=int,
        default=200,
        help="Enforce a length restraint on the the ssODN donor (default: 200), The ssODN donor will be centered on the payload and the recoded region",
        metavar="<int>",
    )
    donor.add_argument(
        "--CheckEnzymes",
        default="",
        help="[dsDNA] Name of Restriction digestion enzymes, separated by '|', to flag and trim, for example BsaI|EcoRI (default: None)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--CustomSeq2Avoid",
        default="",
        help="[dsDNA] Custom sequences, separated by '|', to flag and trim (default: None)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--MinArmLenPostTrim",
        default=0,
        help="[dsDNA] Minimum length of the homology arm after trimming. Set to 0 to turn off trimming (default: 0)",
        type=int,
        metavar="<integer>",
    )
    recoding = parser.add_argument_group('recoding')
    recoding.add_argument(
        "--recoding_off",
        default=False,
        action="store_true",
        help="Turn off *all* recoding",
    )
    recoding.add_argument(
        "--recoding_stop_recut_only",
        default=False,
        action="store_true",
        help="Recode the gRNA recognition site to prevent recut",
    )
    recoding.add_argument(
        "--recoding_full",
        default=False, # if user omit this argument, it will be set to True later in the code near line 317
        action="store_true",
        help="Use full recoding: recode both the gRNA recognition site and the cut-to-insert region (default: on)",
    )
    recoding.add_argument(
        "--recoding_coding_region_only",
        default=False,
        action="store_true",
        help="Only recode the coding region",
    )
    recoding.add_argument(
        "--cfdThres",
        default=0.03,
        help="Threshold that protoSpaceJAM will attempt to lower the recut potential (measured by the CFD score) to (default: 0.03)",
        metavar="<float>",
    )
    recoding.add_argument(
        "--recode_order",
        default="PAM_first",
        help="Prioritize recoding in the PAM or in protospacer, possible values: protospacer_first, PAM_first (default: PAM_first)",
        metavar="<string>",
    )
    misc = parser.add_argument_group('misc.')
    misc.add_argument(
        "--test_mode",
        default=False,
        help="used by the unit tests, not user-oriented",
        metavar="<boolean>",
    )
    misc.add_argument(
        "--guides_only",
        default=False,
        action="store_true",
        help="Only fetch and write guide candidates (skip HDR/recoding/donor design and all genome annotation loads). "
             "Requires coordinate-based input (Chromosome + Coordinate).",
    )
    config = parser.parse_args()
    return config, parser


def main(custom_args=None):
    """
    main function
    custom_args: a dict of arguments to override the default arguments
    """
    try:
        # set up working directory
        if not os.path.exists(os.path.join("genome_files")):
            if os.path.exists(os.path.join("protoSpaceJAM", "genome_files")):
                os.chdir("protoSpaceJAM")

        logging.setLoggerClass(ColoredLogger)
        # logging.basicConfig()
        log = logging.getLogger("protoSpaceJAM")
        log.propagate = False
        log.setLevel(logging.INFO)  # set the level of warning displayed
        # log.setLevel(logging.DEBUG) #set the level of warning displayed

        # configs
        config = vars(parse_args()[0])
        parser = parse_args()[1]

        # apply custom args
        if not custom_args is None and len(custom_args) > 0:
            for c_arg in custom_args:
                config[c_arg] = custom_args[c_arg]

        # Resolve optional payload inputs from files before payload parsing.
        resolve_payload_inputs_from_files(config)
        validate_payload_inputs(config, log)

        # Exit if no arguments provided and not in test mode
        if len(sys.argv)==1 and config["test_mode"] == False:
            print("[Message] Pleases provide the following arguments: --path2csv --outdir")
            print("[Message] To run a quick example: protoSpaceJAM --path2csv input/test_input.csv --outdir output/test\n")

            parser.print_help(sys.stderr)
            sys.exit(1)

        gRNA_num_out = config["num_gRNA_per_design"]
        max_cut2ins_dist = int(config["max_cut2ins_dist"])
        HDR_arm_len = config["HA_len"]
        ssODN_max_size = config["ssODN_max_size"]
        spec_score_flavor = "guideMITScore"
        outdir = config["outdir"]
        reg_penalty = not config["no_regulatory_penalty"]
        guide_source = str(config["guide_source"]).lower()
        chopchop_config = {
            "cmd_template": config["chopchop_cmd_template"],
            "window_padding": config["chopchop_window_padding"],
            "ensembl_timeout": config["ensembl_timeout"],
            "web_base_url": config["chopchop_web_base_url"],
            "web_timeout": config["chopchop_web_timeout"],
            "web_poll_interval": config["chopchop_web_poll_interval"],
            "web_payload_json": config["chopchop_web_payload_json"],
            "debug_dump": config["chopchop_debug_dump"],
            "debug_outdir": outdir,
            "psj_rank_chopchop": config["psj_rank_chopchop"],
            "guides_only": config["guides_only"],
            "scorer_config": {
                "backend": config["specificity_backend"],
                "cmd_template": config["crispor_cmd_template"],
                "max_guides": config["num_gRNA_per_design"],
                "chopchop_proxy_k": config["chopchop_proxy_k"],
                "chopchop_proxy_n": config["chopchop_proxy_n"],
            },
        }
        syn_check_args = {
            "check_enzymes": config["CheckEnzymes"],
            "CustomSeq2Avoid": config["CustomSeq2Avoid"],
            "MinArmLenPostTrim": config["MinArmLenPostTrim"],
        }  # dictionary for multiple synthesis check arguments

        # check recoding args
        assert (
            config["recode_order"] == "protospacer_first"
            or config["recode_order"] == "PAM_first"
        )
        if config["recoding_full"] and any(
            [config["recoding_off"], config["recoding_stop_recut_only"]]
        ):
            sys.exit(
                f"Found conflicts in recoding arguments: --recoding_full cannot be used with --recoding_off or --recoding_stop_recut_only\nplease correct the issue and try again"
            )
        if config["recoding_off"] and config["recoding_stop_recut_only"]:
            sys.exit(
                f"Found conflicts in recoding arguments: --recoding_off cannot be used with --recoding_stop_recut_only\nplease correct the issue and try again"
            )

        # process recoding args
        if (not config["recoding_off"]) and (not config["recoding_stop_recut_only"]) and (not config["recoding_full"]): # if user omit the recoding arguments, it will be set to defaults
            # defaults for SNP mode
            if config["payload_type"] == "SNP":
                config["recoding_stop_recut_only"] = True
                log.info("using default recoding intensity (for SNP payload type): stop recut only")
            else:
                config["recoding_full"] = True
                log.info("using default recoding intensity (for insertion payload type): full, which will recode to stop recut and recoding the cut-to-insert region")
        
        if config["recoding_off"] or config["recoding_stop_recut_only"]:
            config["recoding_full"] = False
            log.info("recoding is turned off by user") if config["recoding_off"] else log.info("recoding intensity is set by the user to prevent recut")

        recoding_args = {
            "recoding_off": config["recoding_off"],
            "recoding_stop_recut_only": config["recoding_stop_recut_only"],
            "recoding_full": config["recoding_full"],
            "cfdThres": float(config["cfdThres"]),
            "recode_order": config["recode_order"],
            "recoding_coding_region_only": config["recoding_coding_region_only"],
        }

        # check donor args
        if not config["Donor_type"] in ["ssODN", "dsDNA"]:
            sys.exit(
                "Donor_type must be ssODN or dsDNA, offending value:"
                + config["Donor_type"]
                + ", please correct the issue and try again"
            )
        if not config["Strand_choice"] in [
            "auto",
            "TargetStrand",
            "NonTargetStrand",
            "CodingStrand",
            "NonCodingStrand",
        ]:
            sys.exit(
                "Strand_choice must be auto,TargetStrand,NonTargetStrand,CodingStrand or NonCodingStrand, offending value:"
                + config["Strand_choice"]
                + ", please correct the issue and try again"
            )

        # check pam
        if not config["pam"].upper() in ["NGG", "NGA", "TTTV"]:
            sys.exit("PAM must be NGG, NGA or TTTV, please correct the issue and try again")
        if guide_source not in ["precomputed", "chopchop"]:
            sys.exit("guide_source must be precomputed or chopchop")
        specificity_backend = str(config["specificity_backend"]).lower()
        if specificity_backend not in ["default", "chopchop_proxy", "crispor"]:
            sys.exit("specificity_backend must be default, chopchop_proxy, or crispor")
        if specificity_backend == "crispor" and str(config["crispor_cmd_template"]).strip() == "":
            sys.exit("specificity_backend=crispor requires --crispor_cmd_template")
        if str(config["input_mode"]).lower() not in ["csv", "direct"]:
            sys.exit("--input_mode must be csv or direct")
        annotation_source = str(config.get("annotation_source", "auto")).lower().strip()
        if annotation_source not in ["auto", "local", "ensembl"]:
            sys.exit("--annotation_source must be auto, local, or ensembl")

        # parse payload
        if config["payload_type"] == "insertion":
            log.info("payload_type is insertion")
            Linker = config["Linker"]
            Tag = config["Tag"]
            if config["payload"] == "":  # no payload override
                if config["Npayload"] == "":  # no Npayload override
                    config["Npayload"] = Tag + Linker
                if config["Cpayload"] == "":  # no Cpayload override
                    config["Cpayload"] = Linker + Tag
            else:  # payload override
                config["Npayload"] = config["payload"]
                config["Cpayload"] = config["payload"]
                config["POSpayload"] = config["payload"]
        elif config["payload_type"] == "SNP":
            log.info("payload_type is SNP")
            if config["SNPpayload"] == "":
                sys.exit("--SNPpayload must be defined for SNP payload type")
        else:
            sys.exit("payload_type must be insertion or SNP, please correct the issue and try again")

        # check if HA_len is too short to satisfy ssODN_max_size
        if ssODN_max_size is not None and config["Donor_type"] == "ssODN":
            max_payload_size = max([len(config["Npayload"]), len(config["Cpayload"])])
            derived_HDR_arm_len = round((ssODN_max_size - max_payload_size) / 2)
            if derived_HDR_arm_len >= HDR_arm_len:
                print(
                    f"HA_len={HDR_arm_len} is to short to meet the requirement of ssODN_max_size={ssODN_max_size}, payload size={max_payload_size}\n ssODN_max_size={ssODN_max_size} requires HA_len = ssODN_max_size- max_payload_size / 2 = {derived_HDR_arm_len}"
                )
                HDR_arm_len = derived_HDR_arm_len + 100
                print(f"HA_len is adjusted to {HDR_arm_len}")

        # check alpha values
        if config["alpha1"] < 0 or config["alpha2"] < 0 or config["alpha3"] < 0:
            sys.exit("alpha values must >= 0, please correct the issue and try again")
        if config["alpha1"] >1 or config["alpha2"] >1  or config["alpha3"] >1 :
            sys.exit("alpha values must be <= 1, please correct the issue and try again")
        if config["alpha1"] == 0 and config["alpha2"] == 0 and config["alpha3"] == 0:
            sys.exit("At least one alpha value must be > 0, please correct the issue and try again")
        alphas = [config["alpha1"], config["alpha2"], config["alpha3"]] 

        # TODO: fix  the potential infinite loop
        # check memory requirement
        enough_mem = test_memory(4200)
        while not enough_mem:  # test if at least 4.2 GB memory is available
            time.sleep(5)  # retry in 5 seconds
            enough_mem = test_memory(4200)

        starttime = datetime.datetime.now()
        freq_dict = dict()
        use_ensembl_annotation_fallback = False
        ensembl_cache_dir = str(config.get("ensembl_cache_dir", "")).strip()
        if ensembl_cache_dir == "":
            ensembl_cache_dir = os.path.join(outdir, "ensembl_cache")
        if annotation_source == "ensembl" or config.get("cache_ensembl_annotation", False):
            mkdir(ensembl_cache_dir)

        loc2file_index = None
        if guide_source == "precomputed":
            # load gRNA info index (mapping of chromosomal location to file parts)
            log.info("loading precomputed gRNA **index to file parts**")
            loc2file_index = read_pickle_files(
                os.path.join(
                    "precomputed_gRNAs",
                    "gRNAs_" + config["pam"].upper(),
                    "gRNA_" + config["genome_ver"],
                    "gRNA.tab.gz.split.BwaMapped.scored",
                    "loc2file_index.pickle",
                )
            )

            elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
            log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")
        else:
            log.info("guide_source=chopchop, skipping precomputed gRNA index loading")

        if config["guides_only"]:
            mkdir(outdir)
            df = load_input_dataframe(config=config, outdir=outdir, log=log)
            guides_out = open(f"{outdir}/guides_from_{guide_source}.csv", "w")
            guides_out.write(
                "Entry,ID,terminus,rank,chopchop_rank,chr,insert_pos,seq,pam,start,end,strand,cut_pos,target_region_label,target_region_start,target_region_end,Eff_scores,MM0,MM1,MM2,MM3,Cut2Ins_dist\n"
            )
            chopchop_raw_all = []
            enst_insert_cache = {}
            ensembl_bundle_cache = {}

            for idx, row in df.iterrows():
                enst_id = str(row.get("Ensembl_ID", "")).strip()
                entry = str(row.get("Entry", idx + 1)).strip()
                chrom = str(row.get("Chromosome", "")).strip()
                coordinate = str(row.get("Coordinate", "")).strip()
                gene_name = str(row.get("Gene_Name", "")).strip()
                this_chopchop_cfg = dict(chopchop_config)
                if gene_name != "":
                    this_chopchop_cfg["gene_input"] = gene_name

                has_coord = (chrom != "" and coordinate.isdigit())
                if (not has_coord) and gene_name == "":
                    log.warning(
                        f"guides_only: skipping entry {entry} ({enst_id}) because Chromosome/Coordinate is missing and Gene_Name is not provided"
                    )
                    continue
                loc_chrom = chrom if chrom != "" else "NA"
                loc_pos = int(coordinate) if has_coord else 1

                try:
                    if guide_source == "chopchop":
                        raw_df = get_chopchop_raw_results(
                            loc=[loc_chrom, loc_pos, 1],
                            dist=max_cut2ins_dist,
                            genome_ver=config["genome_ver"],
                            pam=config["pam"],
                            chopchop_config=this_chopchop_cfg,
                        )
                        raw_df.insert(0, "ID", enst_id if enst_id != "" else gene_name)
                        raw_df.insert(0, "Entry", entry)
                        chopchop_raw_all.append(raw_df)
                        guides_df = convert_chopchop_raw_to_psj(
                            df_raw=raw_df.drop(columns=["Entry", "ID"], errors="ignore"),
                            pam=config["pam"],
                            window_start=1,
                            desired_insert_pos=(int(coordinate) if has_coord else None),
                            default_chr=(chrom if has_coord else None),
                            genome_ver=config["genome_ver"],
                            scorer_config=this_chopchop_cfg.get("scorer_config"),
                        )
                        out_id = enst_id if enst_id != "" else gene_name
                        row_term = str(row.get("Target_terminus", "ALL")).strip().upper()
                        if row_term not in ["N", "C", "ALL", "-"]:
                            row_term = "ALL"

                        if has_coord:
                            write_guide_table_rows(
                                handle=guides_out,
                                entry=entry,
                                enst_id=out_id,
                                terminus="-",
                                ranked_df=guides_df,
                            )
                        else:
                            insert_info = get_insert_positions_from_local_enst(
                                genome_ver=config["genome_ver"],
                                enst_id=enst_id,
                                cache=enst_insert_cache,
                            )
                            if insert_info is None and annotation_source in ["auto", "ensembl"]:
                                insert_info = get_insert_positions_from_ensembl_bundle(
                                    genome_ver=config["genome_ver"],
                                    enst_id=enst_id,
                                    cache=ensembl_bundle_cache,
                                    cache_dir=ensembl_cache_dir,
                                    write_cache=config.get("cache_ensembl_annotation", False),
                                )
                            if insert_info is None:
                                log.warning(
                                    f"guides_only: could not resolve insert positions for {enst_id}; writing with empty insert_pos"
                                )
                                write_guide_table_rows(
                                    handle=guides_out,
                                    entry=entry,
                                    enst_id=out_id,
                                    terminus="-",
                                    ranked_df=guides_df,
                                )
                            else:
                                n_ins, c_ins, ann_chr = insert_info
                                if ann_chr and ("chr" in guides_df.columns):
                                    guides_df["chr"] = ann_chr
                                if row_term in ["N", "ALL"]:
                                    df_n = guides_df.copy()
                                    df_n["Insert_pos"] = n_ins
                                    write_guide_table_rows(
                                        handle=guides_out,
                                        entry=entry,
                                        enst_id=out_id,
                                        terminus="N",
                                        ranked_df=df_n,
                                    )
                                if row_term in ["C", "ALL"]:
                                    df_c = guides_df.copy()
                                    df_c["Insert_pos"] = c_ins
                                    write_guide_table_rows(
                                        handle=guides_out,
                                        entry=entry,
                                        enst_id=out_id,
                                        terminus="C",
                                        ranked_df=df_c,
                                    )
                    elif has_coord:
                        guides_df = get_gRNAs_near_loc(
                            loc=[chrom, int(coordinate), 1],
                            dist=max_cut2ins_dist,
                            loc2file_index=loc2file_index,
                            genome_ver=config["genome_ver"],
                            pam=config["pam"],
                            guide_source=guide_source,
                            chopchop_config=this_chopchop_cfg,
                        )
                        if str(config["specificity_backend"]).lower() != "default" and "guideMITScore" in guides_df.columns:
                            guides_df = guides_df.sort_values(
                                "guideMITScore", ascending=False
                            ).reset_index(drop=True)
                        write_guide_table_rows(
                            handle=guides_out,
                            entry=entry,
                            enst_id=enst_id if enst_id != "" else gene_name,
                            terminus="-",
                            ranked_df=guides_df,
                        )
                except Exception as e:
                    log.warning(f"guides_only: failed entry {entry} ({enst_id}): {e}")

            guides_out.close()
            if guide_source == "chopchop" and len(chopchop_raw_all) > 0:
                raw_out_path = os.path.join(outdir, "guides_from_chopchop_web_fields.csv")
                pd.concat(chopchop_raw_all, ignore_index=True).to_csv(raw_out_path, index=False)
                log.info(f"guides_only raw CHOPCHOP fields written to {raw_out_path}")
            log.info(
                f"guides_only completed. Guide table written to {os.path.join(outdir, f'guides_from_{guide_source}.csv')}"
            )
            return

        loc2posType_path = os.path.join(
            "genome_files",
            "parsed_gff3",
            config["genome_ver"],
            "loc2posType.pickle",
        )
        enst_info_idx_path = os.path.join(
            "genome_files",
            "parsed_gff3",
            config["genome_ver"],
            "ENST_info",
            "ENST_info_index.pickle",
        )
        codon_idx_path = os.path.join(
            "genome_files",
            "parsed_gff3",
            config["genome_ver"],
            "ENST_codonPhases",
            "ENST_codonPhase_index.pickle",
        )
        local_annotations_available = (
            os.path.isfile(loc2posType_path)
            and os.path.isfile(enst_info_idx_path)
            and os.path.isfile(codon_idx_path)
        )
        if annotation_source == "local" and (not local_annotations_available):
            sys.exit(
                "annotation_source=local but required local annotation pickles are missing "
                f"(expected: {loc2posType_path}, {enst_info_idx_path}, {codon_idx_path})"
            )
        if annotation_source == "ensembl":
            use_ensembl_annotation_fallback = True
        elif annotation_source == "auto":
            use_ensembl_annotation_fallback = (not local_annotations_available)

        if use_ensembl_annotation_fallback:
            loc2posType = {}
            log.warning(
                "using Ensembl on-the-fly transcript annotation source (optional local cache enabled)"
            )
        else:
            # load chr location to type (e.g. UTR, cds, exon/intron junction) mappings
            log.info(
                "loading the mapping of chromosomal location to type (e.g. UTR, cds, exon/intron junction)"
            )
            loc2posType = read_pickle_files(loc2posType_path)
            elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
            log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # load gene model info
        # log.info("loading gene model info")
        # ENST_info = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_info.pickle"))

        if use_ensembl_annotation_fallback:
            ENST_info_index = {}
        else:
            # load the mapping of ENST to the info file part
            log.info("loading gene model info **index to file parts**")
            ENST_info_index = read_pickle_files(
                os.path.join(
                    "genome_files",
                    "parsed_gff3",
                    config["genome_ver"],
                    "ENST_info",
                    "ENST_info_index.pickle",
                )
            )
            elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
            log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # load codon phase index
        # log.info("loading codon phase info")
        # ENST_PhaseInCodon = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhase.pickle"))

        if use_ensembl_annotation_fallback:
            ENST_PhaseInCodon_index = {}
        else:
            log.info("loading codon phase info **index to file parts**")
            ENST_PhaseInCodon_index = read_pickle_files(
                os.path.join(
                    "genome_files",
                    "parsed_gff3",
                    config["genome_ver"],
                    "ENST_codonPhases",
                    "ENST_codonPhase_index.pickle",
                )
            )

        # report time used
        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # report the number of ENSTs which has ATG at the end of the exon
        # ExonEnd_ATG_count,ExonEnd_ATG_list = count_ATG_at_exonEnd(ENST_info)

        # open log files
        mkdir(outdir)
        mkdir(os.path.join(outdir, "genbank_files"))
        # Legacy CFD and fiveUTR side outputs are kept in-memory to avoid cluttering the output directory.
        recut_CFD_all = io.StringIO()
        recut_CFD_fail = io.StringIO()
        csvout_N = io.StringIO()
        csvout_C = io.StringIO()
        csvout_header = "ID,cfd1,cfd2,cfd3,cfd4,cfdScan,cfdScanNoRecode,cfd_max\n"
        csvout_N.write(csvout_header)
        csvout_C.write(csvout_header)
        fiveUTR_log = io.StringIO()
        ha_out = open(os.path.join(outdir, "homology_arms.csv"), "w")
        ha_out.write(
            "Entry,ID,terminus,design_rank,gRNA_name,gRNA_seq,insert_pos,left_HA,payload,right_HA,donor_final\n"
        )

        # open result file and write header
        csvout_res = open(f"{outdir}/result.csv", "w")
        csvout_res.write(
            f"Entry,ID,chr,transcript_type,name,terminus,design_rank,gRNA_name,gRNA_seq,PAM,gRNA_start,gRNA_end,gRNA_cut_pos,edit_pos,distance_between_cut_and_edit(cut_pos-insert_pos),chopchop_rank,target_region_label,target_region_start,target_region_end,cfd_before_recoding,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,name_of_DNA_donor,DNA donor,name_of_trimmed_DNA_Donor,trimmed_DNA_donor,effective_HA_len,synthesis_problems,cutPos2nearestOffLimitJunc,strand(gene/gRNA/donor)\n"
        )   #"Entry,ID,chr,transcript_type,name,terminus,gRNA_seq,PAM,gRNA_start,gRNA_end,gRNA_cut_pos,edit_pos,distance_between_cut_and_edit(cut pos - insert pos),specificity_score,specificity_weight,distance_weight,position_weight,final_weight,cfd_before_recoding,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,DNA donor,effective_HA_len,synthesis_problems,cutPos2nearestOffLimitJunc,strand(gene/gRNA/donor)\n"

        # Legacy GenoPrimer export is kept in-memory for now.
        csvout_res2 = io.StringIO()
        csvout_res2.write(f"Entry,ref,chr,coordinate,ID,geneSymbol\n")
        guides_out = open(f"{outdir}/guides_from_{guide_source}.csv", "w")
        guides_out.write(
            "Entry,ID,terminus,rank,chopchop_rank,chr,insert_pos,seq,pam,start,end,strand,cut_pos,target_region_label,target_region_start,target_region_end,Eff_scores,MM0,MM1,MM2,MM3,Cut2Ins_dist\n"
        )

        # dataframes to store best gRNAs
        best_start_gRNAs = pd.DataFrame()
        best_stop_gRNAs = pd.DataFrame()

        # logging cfd score, failed gRNAs etc
        start_info = info()
        stop_info = info()

        # load ENST list (from CSV) or create a one-row input table from direct CLI args
        df = load_input_dataframe(config=config, outdir=outdir, log=log)
        guides_df = load_guides_dataframe(config=config, log=log)

        # loop through each entry in the input csv file
        transcript_count = 0
        protein_coding_transcripts_count = 0
        target_terminus = "None"
        target_coordinate = "None"
        ENST_in_db = False
        ENST_design_counts = {} #used to keep track of number of designs for each ENST
        ensembl_info_cache = {}
        Entry = 0 # Entry is defined by the portal, and if not using the portal, it is just the index of the row
        total_entries = df.shape[0]
        skipped_count = 0
        for index, row in df.iterrows():
            target_terminus = "ALL"
            chrom = ""
            coordinate = ""
            ENST_ID = row["Ensembl_ID"]
            if isinstance(ENST_ID, str):
                ENST_ID =  ENST_ID.rstrip().lstrip()
            gene_name = ""
            if "Gene_Name" in df.columns and row.isnull().get("Gene_Name", True) == False:
                gene_name = str(row["Gene_Name"]).rstrip().lstrip()
            this_chopchop_cfg = dict(chopchop_config)
            if gene_name != "":
                this_chopchop_cfg["gene_input"] = gene_name
            if "Target_terminus" in df.columns and row.isnull()["Target_terminus"] == False:
                target_terminus = row["Target_terminus"].rstrip().lstrip().upper()

                if (
                    target_terminus != "N"
                    and target_terminus != "C"
                    and target_terminus != "ALL"
                    and target_terminus != ""
                ):
                    sys.exit(f"invalid target terminus: {target_terminus}")

            # determine if the input is ENST-based or coordinate-based
            ENST_based, coordinate_based, coordinate_without_ENST = False, False, False
            if "Chromosome" in df.columns and "Coordinate" in df.columns and row.isnull()["Chromosome"] == False and row.isnull()["Coordinate"] == False:
                chrom = str(row["Chromosome"]).rstrip().lstrip()
                coordinate = str(row["Coordinate"]).rstrip().lstrip()
                if (
                    chrom != ""
                    and coordinate != ""
                    and coordinate.isdigit()
                ):
                    coordinate_based = True
                    # check if ENST is provided,
                    if (not use_ensembl_annotation_fallback) and (not ENST_ID in ENST_info_index.keys()):
                        coordinate_without_ENST = True # support for coordinate-based design without ENST
                        # create a dummy ENST_ID for the coordinate-based design without ENST
                        if len(ENST_PhaseInCodon_index) > 0:
                            ENST_ID = next(iter(ENST_PhaseInCodon_index))
                    if use_ensembl_annotation_fallback and ENST_ID == "":
                        coordinate_without_ENST = True

            # if coordinate_based didn't check out, revert to ENST-based
            if coordinate_based == False:
                ENST_based = True
                if not target_terminus in ["N","C","ALL"]:
                    target_terminus = "ALL"

            if "Entry" in df.columns:
                try:
                    Entry = str(row["Entry"]).rstrip().lstrip()
                    ENST_in_db = True # this bool indicates that the input is from the portal
                except:
                    ENST_in_db = False

            # check if ENST_ID is in the database
            if (   (not use_ensembl_annotation_fallback) and (not ENST_ID in ENST_info_index.keys()) 
                and not coordinate_without_ENST
                ):
                if not ENST_in_db: # case where not using the portal
                    Entry += 1
                log.warning(
                    f"skipping {ENST_ID} b/c transcript is not in the annotated ENST collection (excluding those on chr_patch_hapl_scaff)"
                )
                genome_ver = config["genome_ver"]
                csvout_res.write(
                    f"{Entry},{ENST_ID},ERROR: this ID was not found in the genome {genome_ver}, most likely this ID was deprecated\n"
                )
                continue

            # check if codon phase info exists
            if (   (not use_ensembl_annotation_fallback) and (not ENST_ID in ENST_PhaseInCodon_index.keys()) 
                and not coordinate_without_ENST):
                if not ENST_in_db:
                    Entry += 1
                log.warning(
                    f"skipping {ENST_ID} b/c transcript has no codon phase information"
                )
                csvout_res.write(
                    f"{Entry},{ENST_ID},ERROR: this ID is either not protein-coding and/or has no codon phase information\n"
                )
                continue

            # load the ENST_info for current ID
            if use_ensembl_annotation_fallback:
                enst_norm = _normalize_enst_id(ENST_ID)
                bundle = get_ensembl_annotation_bundle(
                    enst_id=ENST_ID,
                    genome_ver=config["genome_ver"],
                    cache=ensembl_info_cache,
                    cache_dir=ensembl_cache_dir,
                    write_cache=config.get("cache_ensembl_annotation", False),
                )
                if bundle is None or "ENST_info" not in bundle or enst_norm not in bundle["ENST_info"]:
                    log.warning(
                        f"skipping {ENST_ID}: could not retrieve transcript annotation from Ensembl"
                    )
                    csvout_res.write(
                        f"{Entry},{ENST_ID},ERROR: unable to retrieve transcript annotation from Ensembl\n"
                    )
                    continue
                ENST_info = bundle["ENST_info"]
                ENST_ID = enst_norm
                if "loc2posType" in bundle and len(bundle["loc2posType"]) > 0:
                    loc2posType = deepmerge(loc2posType, bundle["loc2posType"])
            else:
                part = ENST_info_index[ENST_ID]
                ENST_info = read_pickle_files(
                    os.path.join(
                        "genome_files",
                        "parsed_gff3",
                        config["genome_ver"],
                        "ENST_info",
                        f"ENST_info_part{part}.pickle",
                    )
                )

            transcript_type = ENST_info[ENST_ID].description.split("|")[1]
            # if transcript_type == "protein_coding": # and ENST_ID == "ENST00000398165":
            # if not ENST_ID in ExonEnd_ATG_list: # only process edge cases in which genes with ATG are at the end of exons
            #     continue
            log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")

            resolved_site = None
            preferred_region = _clean_optional_str(row.get("Preferred_region", ""))
            preferred_exon = _clean_optional_str(row.get("Preferred_exon", ""))
            preferred_anchor = _clean_optional_str(row.get("Preferred_anchor", ""))
            preferred_offset = _clean_optional_str(row.get("Preferred_offset", ""))
            preferred_mode_requested = any(
                val != "" for val in [preferred_region, preferred_exon, preferred_anchor, preferred_offset]
            )
            if (not coordinate_based) and preferred_mode_requested:
                try:
                    resolved_site = resolve_preferred_insert_site(
                        ENST_info=ENST_info,
                        ENST_ID=ENST_ID,
                        preferred_region=preferred_region,
                        preferred_exon=preferred_exon,
                        preferred_anchor=preferred_anchor,
                        preferred_offset=preferred_offset,
                    )
                except Exception as e:
                    csvout_res.write(
                        f"{Entry},{ENST_ID},ERROR: could not resolve Preferred_region: {str(e).replace(',', ';')}\n"
                    )
                    continue
                chrom = str(resolved_site["chrom"])
                coordinate = str(resolved_site["coordinate"])
                coordinate_based = True
                ENST_based = False
                log.info(
                    f"resolved preferred region {resolved_site['region_label']} ({resolved_site['anchor']}, offset={resolved_site['offset']}) to {chrom}:{coordinate}"
                )

            if hasattr(ENST_info[ENST_ID], "name"):
                name = ENST_info[ENST_ID].name
            else:
                name = ""
            row_prefix = f"{ENST_ID},{ENST_info[ENST_ID].chr},{transcript_type},{name}"
            if coordinate_without_ENST:
                row_prefix = f",,,"

            # get codon_phase information for current ENST
            ENST_PhaseInCodon = {}
            if use_ensembl_annotation_fallback:
                if "ENST_PhaseInCodon" in bundle and len(bundle["ENST_PhaseInCodon"]) > 0:
                    ENST_PhaseInCodon = deepmerge(ENST_PhaseInCodon, bundle["ENST_PhaseInCodon"])
            elif ENST_ID in ENST_PhaseInCodon_index:
                file_parts_list = ENST_PhaseInCodon_index[ENST_ID]
                for part in file_parts_list:
                    Codon_phase_dict = read_pickle_files(
                        os.path.join(
                            "genome_files",
                            "parsed_gff3",
                            config["genome_ver"],
                            "ENST_codonPhases",
                            f"ENST_codonPhase_part{str(part)}.pickle",
                        )
                    )
                    ENST_PhaseInCodon = deepmerge(
                        ENST_PhaseInCodon, Codon_phase_dict
                    )  # merge file parts if ENST codon info is split among fileparts

            ######################################
            # best gRNA for a specific coordinate#
            ######################################
            if coordinate_based == True:
                if not ENST_in_db:
                    Entry += 1
                csvout_N.write(ENST_ID)
                #check if the coordinate is in the ENST
                coord_in_ENST = ENST_info[ENST_ID].chr == chrom and min(ENST_info[ENST_ID].span_start, ENST_info[ENST_ID].span_end) <= int(coordinate) <=  max(ENST_info[ENST_ID].span_start, ENST_info[ENST_ID].span_end)

                if not coordinate_without_ENST and not coord_in_ENST: # if the coordinate is not in the ENST (and user provided a ENST_ID to go with the coordinate)
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: provided genomic coordinates are not in the {ENST_ID}\n")
                else:                    
                    # get gRNAs
                    if guides_df is not None:
                        ranked_df_gRNAs_target_pos = select_guides_from_csv(
                            guides_df=guides_df,
                            entry=Entry,
                            enst_id=ENST_ID,
                            terminus="-",
                            chrom=chrom,
                            coordinate=coordinate,
                            max_abs_cut2ins=max_cut2ins_dist,
                        )
                        ranked_df_gRNAs_target_pos = rerank_guides_csv_with_psj(
                            guides_df=ranked_df_gRNAs_target_pos,
                            enst_id=ENST_ID,
                            chrom=chrom,
                            insert_pos=int(coordinate),
                            enst_strand=ENST_info[ENST_ID].features[0].strand,
                            terminus_type="start",
                            loc2posType=loc2posType,
                            spec_score_flavor=spec_score_flavor,
                            reg_penalty=reg_penalty,
                            alphas=alphas,
                            max_abs_cut2ins=max_cut2ins_dist,
                            specificity_backend=config["specificity_backend"],
                            chopchop_proxy_k=config["chopchop_proxy_k"],
                            chopchop_proxy_n=config["chopchop_proxy_n"],
                            preserve_input_order=str(config["specificity_backend"]).lower() == "default",
                        )
                    else:
                        try:
                            ranked_df_gRNAs_target_pos = get_gRNAs_target_coordinate(
                                ENST_ID=ENST_ID,
                                chrom = chrom,
                                pos = int(coordinate),
                                ENST_info=ENST_info,
                                freq_dict=freq_dict,
                                loc2file_index=loc2file_index,
                                loc2posType=loc2posType,
                                dist=max_cut2ins_dist,
                                genome_ver=config["genome_ver"],
                                pam=config["pam"],
                                spec_score_flavor=spec_score_flavor,
                                reg_penalty=reg_penalty,
                                alphas = alphas,
                                guide_source=guide_source,
                                chopchop_config=this_chopchop_cfg,
                            )
                        except Exception as e:
                            csvout_res.write(
                                f"{Entry},{ENST_ID},ERROR: guide retrieval failed: {str(e).replace(',', ';')}\n"
                            )
                            continue
                    region_label = "coordinate"
                    region_start = None
                    region_end = None
                    if resolved_site is not None:
                        region_label = resolved_site.get("region_label", region_label)
                        region_start = resolved_site.get("region_start", None)
                        region_end = resolved_site.get("region_end", None)
                    ranked_df_gRNAs_target_pos = _annotate_guides_with_output_region(
                        ranked_df_gRNAs_target_pos,
                        region_label=region_label,
                        region_start=region_start,
                        region_end=region_end,
                    )
                    write_guide_table_rows(guides_out, Entry, ENST_ID, "-", ranked_df_gRNAs_target_pos)

                    if ranked_df_gRNAs_target_pos.empty == True:
                        csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                    for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_target_pos.shape[0]])):
                        current_gRNA = ranked_df_gRNAs_target_pos.iloc[[i]]

                        # get HDR template
                        try:
                            HDR_template = get_HDR_template(
                                df=current_gRNA,
                                ENST_info=ENST_info,
                                type="start", # this borrowed option specifies that the edit is immediately after the coordinate
                                ENST_PhaseInCodon=ENST_PhaseInCodon,
                                loc2posType=loc2posType,
                                genome_ver=config["genome_ver"],
                                HDR_arm_len=HDR_arm_len,
                                payload_type=config["payload_type"],
                                tag=config["POSpayload"],
                                SNP_payload=config["SNPpayload"],
                                ssODN_max_size=ssODN_max_size,
                                Donor_type=config["Donor_type"],
                                Strand_choice=config["Strand_choice"],
                                recoding_args=recoding_args,
                                syn_check_args=syn_check_args,
                                coordinate_without_ENST = coordinate_without_ENST,
                            )
                        except Exception as e:
                            print("Unexpected error:", str(sys.exc_info()))
                            traceback.print_exc()
                            print("additional information:", e)
                            PrintException()
                            guide_seq = current_gRNA["seq"].values[0] if "seq" in current_gRNA.columns else ""
                            csvout_res.write(
                                f"{Entry},{ENST_ID},ERROR: HDR design failed for guide {guide_seq}: {str(e).replace(',', ';')}\n"
                            )
                            continue

                        # append the best gRNA to the final df
                        if i == 0:
                            best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                        # append cfd score to list for plotting
                        pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                        cfd1 = ""
                        if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                            cfd1 = HDR_template.cfd_score_post_mut_ins
                        if not hasattr(HDR_template, "cfd_score_post_mut2"):
                            cfd2 = cfd1
                        else:
                            cfd2 = HDR_template.cfd_score_post_mut2
                        if not hasattr(HDR_template, "cfd_score_post_mut3"):
                            cfd3 = cfd2
                        else:
                            cfd3 = HDR_template.cfd_score_post_mut3
                        if not hasattr(HDR_template, "cfd_score_post_mut4"):
                            cfd4 = cfd3
                        else:
                            cfd4 = HDR_template.cfd_score_post_mut4
                        cfd_scan = 0
                        cfd_scan_no_recode = 0
                        if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                            cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                            cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                        cfdfinal = HDR_template.final_cfd

                        strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                        # write csv
                        (
                            spec_score,
                            seq,
                            pam,
                            s,
                            e,
                            cut2ins_dist,
                            spec_weight,
                            dist_weight,
                            pos_weight,
                            final_weight,
                        ) = get_res(current_gRNA, spec_score_flavor)
                        chopchop_rank = current_gRNA["chopchop_rank"].values[0] if "chopchop_rank" in current_gRNA.columns else ""
                        target_region_label = current_gRNA["target_region_label"].values[0] if "target_region_label" in current_gRNA.columns else ""
                        target_region_start = current_gRNA["target_region_start"].values[0] if "target_region_start" in current_gRNA.columns else ""
                        target_region_end = current_gRNA["target_region_end"].values[0] if "target_region_end" in current_gRNA.columns else ""
                        cut_in_target_region = current_gRNA["cut_in_target_region"].values[0] if "cut_in_target_region" in current_gRNA.columns else ""
                        donor = HDR_template.Donor_final
                        donor_trimmed = "N/A for ssODN"
                        if config["Donor_type"] == "dsDNA":
                            donor_trimmed = HDR_template.Donor_final
                            donor = HDR_template.Donor_pretrim
                            
                        # gRNA and donor names
                        # Include gRNA rank in the name when multiple gRNAs are requested
                        gRNA_rank_suffix = f"_rank{i+1}" if gRNA_num_out > 1 else ""

                        if coordinate_without_ENST:
                            gRNA_name = f"{ENST_ID}_coord_gRNA_entry{Entry}{gRNA_rank_suffix}"
                            donor_name = f"{ENST_ID}_coord_donor_entry{Entry}{gRNA_rank_suffix}"
                        else:
                            chrom_coord = f"{chrom}_{coordinate}"
                            ENST_design_counts[chrom_coord] = ENST_design_counts.get(chrom_coord, 0) + 1
                            gRNA_name = f"{chrom_coord}_gRNA_entry{Entry}{gRNA_rank_suffix}"
                            donor_name = f"{chrom_coord}_donor_entry{Entry}{gRNA_rank_suffix}"

                        if config["Donor_type"] == "dsDNA":
                            # replace donor_ with donor_trimmed_
                            donor_trimmed_name = donor_name.replace("donor_", "donor_trimmed_")
                            #donor_trimmed_name = f"{ENST_ID}_donor_trimmed_design{ENST_design_counts[ENST_ID]}"
                        else:
                            donor_trimmed_name = "N/A for ssODN"

                        gRNA_cut_pos = (
                            HDR_template.CutPos
                        )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                        insert_pos = HDR_template.InsPos
                        if config["recoding_off"]:
                            csvout_res.write(
                                f"{Entry},{row_prefix},-,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                            )
                            csvout_res2.write( f"{Entry},"+
                                config["genome_ver"]
                                + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                            )
                        else:
                            if not isinstance(cfd4, float):
                                cfd4 = ""
                            csvout_res.write(
                                f"{Entry},{row_prefix},-,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                            )
                            csvout_res2.write( f"{Entry},"+
                                config["genome_ver"]
                                + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                            )

                        # donor features
                        donor_features = HDR_template.Donor_features
                        write_ha_csv_row(ha_out, Entry, ENST_ID, "-", i+1, gRNA_name, seq, insert_pos, HDR_template)
                        # write genbank file
                        with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                            write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])
                        # write anothergenbank file without payload
                        with open(os.path.join(outdir, "genbank_files", f"{donor_name}_gRNAonly_noPayload.gb"), "w") as gb_handle:
                            write_genbank_gRNAonly_noPayload(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])

                        # write log
                        this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n    donor before any recoding:{HDR_template.Donor_vanillia}\n     donor after all recoding:{HDR_template.Donor_postMut}\ndonor centered(if applicable):{HDR_template.Donor_final}\n          donor (best strand):{HDR_template.Donor_final}\n\n"
                        #this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                        recut_CFD_all.write(this_log)
                        if HDR_template.final_cfd > 0.03:
                            recut_CFD_fail.write(this_log)

                        if hasattr(HDR_template, "info_phase4_5UTR"):
                            fiveUTR_log.write(
                                f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                            )
                        if hasattr(HDR_template, "info_phase5_5UTR"):
                            fiveUTR_log.write(
                                f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                            )


            ###################################
            # best start gRNA and HDR template#
            ###################################
            if ENST_based == True and (target_terminus == "ALL" or target_terminus == "N"):
                if not ENST_in_db:
                    Entry += 1
                csvout_N.write(ENST_ID)
                n_insert_pos = None
                c_insert_pos = None
                try:
                    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
                    n_insert_pos = int(get_end_pos_of_ATG(ATG_loc)[1])
                    c_insert_pos = int(get_start_pos_of_stop(stop_loc)[1])
                except Exception:
                    n_insert_pos = None
                    c_insert_pos = None
                # get gRNAs
                if guides_df is not None:
                    ranked_df_gRNAs_ATG = select_guides_from_csv(
                        guides_df=guides_df,
                        entry=Entry,
                        enst_id=ENST_ID,
                        terminus="N",
                        fallback_insert_pos=n_insert_pos,
                        fallback_chr=ENST_info[ENST_ID].chr if ENST_ID in ENST_info else None,
                        max_abs_cut2ins=max_cut2ins_dist,
                    )
                    ranked_df_gRNAs_ATG = rerank_guides_csv_with_psj(
                        guides_df=ranked_df_gRNAs_ATG,
                        enst_id=ENST_ID,
                        chrom=ENST_info[ENST_ID].chr,
                        insert_pos=n_insert_pos if n_insert_pos is not None else 0,
                        enst_strand=ENST_info[ENST_ID].features[0].strand,
                        terminus_type="start",
                        loc2posType=loc2posType,
                        spec_score_flavor=spec_score_flavor,
                        reg_penalty=reg_penalty,
                        alphas=alphas,
                        max_abs_cut2ins=max_cut2ins_dist,
                        specificity_backend=config["specificity_backend"],
                        chopchop_proxy_k=config["chopchop_proxy_k"],
                        chopchop_proxy_n=config["chopchop_proxy_n"],
                        preserve_input_order=str(config["specificity_backend"]).lower() == "default",
                    )
                    ranked_df_gRNAs_stop = pd.DataFrame()
                else:
                    try:
                        ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(
                            ENST_ID=ENST_ID,
                            ENST_info=ENST_info,
                            freq_dict=freq_dict,
                            loc2file_index=loc2file_index,
                            loc2posType=loc2posType,
                            dist=max_cut2ins_dist,
                            genome_ver=config["genome_ver"],
                            pam=config["pam"],
                            spec_score_flavor=spec_score_flavor,
                            reg_penalty=reg_penalty,
                            alphas = alphas,
                            guide_source=guide_source,
                            chopchop_config=this_chopchop_cfg,
                        )
                    except Exception as e:
                        csvout_res.write(
                            f"{Entry},{ENST_ID},ERROR: guide retrieval failed: {str(e).replace(',', ';')}\n"
                        )
                        continue
                write_guide_table_rows(guides_out, Entry, ENST_ID, "N", ranked_df_gRNAs_ATG)

                if ranked_df_gRNAs_ATG.empty == True:
                    start_info.failed.append(ENST_ID)
                    csvout_N.write(",,,,,\n")
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_ATG.shape[0]])):
                    # if best_start_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_start_gRNA = best_start_gRNA[best_start_gRNA["CSS"] == best_start_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_start_gRNA = best_start_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_ATG.iloc[[i]]

                    # get HDR template
                    try:
                        HDR_template = get_HDR_template(
                            df=current_gRNA,
                            ENST_info=ENST_info,
                            type="start",
                            ENST_PhaseInCodon=ENST_PhaseInCodon,
                            loc2posType=loc2posType,
                            genome_ver=config["genome_ver"],
                            HDR_arm_len=HDR_arm_len,
                            payload_type=config["payload_type"],
                            tag=config["Npayload"],
                            SNP_payload=config["SNPpayload"],
                            ssODN_max_size=ssODN_max_size,
                            Donor_type=config["Donor_type"],
                            Strand_choice=config["Strand_choice"],
                            recoding_args=recoding_args,
                            syn_check_args=syn_check_args,
                            coordinate_without_ENST = coordinate_without_ENST,
                        )
                    except Exception as e:
                        print("Unexpected error:", str(sys.exc_info()))
                        traceback.print_exc()
                        print("additional information:", e)
                        PrintException()
                        guide_seq = current_gRNA["seq"].values[0] if "seq" in current_gRNA.columns else ""
                        csvout_res.write(
                            f"{Entry},{ENST_ID},ERROR: HDR design failed for guide {guide_seq}: {str(e).replace(',', ';')}\n"
                        )
                        continue

                    # append the best gRNA to the final df
                    if i == 0:
                        best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                    # append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                        cfd1 = HDR_template.cfd_score_post_mut_ins
                    if not hasattr(HDR_template, "cfd_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cfd_score_post_mut2
                    if not hasattr(HDR_template, "cfd_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cfd_score_post_mut3
                    if not hasattr(HDR_template, "cfd_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cfd_score_post_mut4
                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                    # write csv
                    (
                        spec_score,
                        seq,
                        pam,
                        s,
                        e,
                        cut2ins_dist,
                        spec_weight,
                        dist_weight,
                        pos_weight,
                        final_weight,
                    ) = get_res(current_gRNA, spec_score_flavor)

                    donor = HDR_template.Donor_final
                    donor_trimmed = "N/A for ssODN"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed = HDR_template.Donor_final
                        donor = HDR_template.Donor_pretrim

                    # gRNA and donor names
                    # Include gRNA rank in the name when multiple gRNAs are requested
                    gRNA_rank_suffix = f"_rank{i+1}" if gRNA_num_out > 1 else ""

                    ENST_design_counts[ENST_ID] = ENST_design_counts.get(ENST_ID, 0) + 1
                    gRNA_name = f"{ENST_ID}_N_gRNA_entry{Entry}{gRNA_rank_suffix}"
                    donor_name = f"{ENST_ID}_N_donor_entry{Entry}{gRNA_rank_suffix}"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed_name = donor_name.replace("donor_", "donor_trimmed_")
                        #donor_trimmed_name = f"{ENST_ID}_donor_trimmed_{ENST_design_counts[ENST_ID]}"
                    else:
                        donor_trimmed_name = "N/A for ssODN"
                            
                    gRNA_cut_pos = (
                        HDR_template.CutPos
                    )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_N.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        csvout_res.write(
                            f"{Entry},{row_prefix},N,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )
                    else:
                        csvout_N.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        if not isinstance(cfd4, float):
                            cfd4 = ""
                        csvout_res.write(
                            f"{Entry},{row_prefix},N,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )

                    # donor features
                    donor_features = HDR_template.Donor_features
                    write_ha_csv_row(ha_out, Entry, ENST_ID, "N", i+1, gRNA_name, seq, insert_pos, HDR_template)
                    # write genbank file
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                        write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])
                    # write anothergenbank file without payload
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}_gRNAonly_noPayload.gb"), "w") as gb_handle:
                        write_genbank_gRNAonly_noPayload(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])

                    # write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n    donor before any recoding:{HDR_template.Donor_vanillia}\n     donor after all recoding:{HDR_template.Donor_postMut}\ndonor centered(if applicable):{HDR_template.Donor_final}\n          donor (best strand):{HDR_template.Donor_final}\n\n"
                    #this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template, "info_phase4_5UTR"):
                        fiveUTR_log.write(
                            f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                        )
                    if hasattr(HDR_template, "info_phase5_5UTR"):
                        fiveUTR_log.write(
                            f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                        )

            ##################################
            # best stop gRNA and HDR template#
            ##################################
            if ENST_based == True and (target_terminus == "ALL" or target_terminus == "C"):
                if not ENST_in_db:
                    Entry += 1
                csvout_C.write(ENST_ID)
                n_insert_pos = None
                c_insert_pos = None
                try:
                    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
                    n_insert_pos = int(get_end_pos_of_ATG(ATG_loc)[1])
                    c_insert_pos = int(get_start_pos_of_stop(stop_loc)[1])
                except Exception:
                    n_insert_pos = None
                    c_insert_pos = None

                # get gRNAs
                if guides_df is not None:
                    ranked_df_gRNAs_ATG = pd.DataFrame()
                    ranked_df_gRNAs_stop = select_guides_from_csv(
                        guides_df=guides_df,
                        entry=Entry,
                        enst_id=ENST_ID,
                        terminus="C",
                        fallback_insert_pos=c_insert_pos,
                        fallback_chr=ENST_info[ENST_ID].chr if ENST_ID in ENST_info else None,
                        max_abs_cut2ins=max_cut2ins_dist,
                    )
                    ranked_df_gRNAs_stop = rerank_guides_csv_with_psj(
                        guides_df=ranked_df_gRNAs_stop,
                        enst_id=ENST_ID,
                        chrom=ENST_info[ENST_ID].chr,
                        insert_pos=c_insert_pos if c_insert_pos is not None else 0,
                        enst_strand=ENST_info[ENST_ID].features[0].strand,
                        terminus_type="stop",
                        loc2posType=loc2posType,
                        spec_score_flavor=spec_score_flavor,
                        reg_penalty=reg_penalty,
                        alphas=alphas,
                        max_abs_cut2ins=max_cut2ins_dist,
                        specificity_backend=config["specificity_backend"],
                        chopchop_proxy_k=config["chopchop_proxy_k"],
                        chopchop_proxy_n=config["chopchop_proxy_n"],
                        preserve_input_order=str(config["specificity_backend"]).lower() == "default",
                    )
                else:
                    try:
                        ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(
                            ENST_ID=ENST_ID,
                            ENST_info=ENST_info,
                            freq_dict=freq_dict,
                            loc2file_index=loc2file_index,
                            loc2posType=loc2posType,
                            dist=max_cut2ins_dist,
                            genome_ver=config["genome_ver"],
                            pam=config["pam"],
                            spec_score_flavor=spec_score_flavor,
                            reg_penalty=reg_penalty,
                            alphas = alphas,
                            guide_source=guide_source,
                            chopchop_config=this_chopchop_cfg,
                        )
                    except Exception as e:
                        csvout_res.write(
                            f"{Entry},{ENST_ID},ERROR: guide retrieval failed: {str(e).replace(',', ';')}\n"
                        )
                        continue
                write_guide_table_rows(guides_out, Entry, ENST_ID, "C", ranked_df_gRNAs_stop)

                if ranked_df_gRNAs_stop.empty == True:
                    stop_info.failed.append(ENST_ID)
                    csvout_C.write(",,,,,\n")
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_stop.shape[0]])):
                    # if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_stop.iloc[[i]]

                    # get HDR template
                    try:
                        HDR_template = get_HDR_template(
                            df=current_gRNA,
                            ENST_info=ENST_info,
                            type="stop",
                            ENST_PhaseInCodon=ENST_PhaseInCodon,
                            loc2posType=loc2posType,
                            HDR_arm_len=HDR_arm_len,
                            genome_ver=config["genome_ver"],
                            payload_type=config["payload_type"],
                            tag=config["Cpayload"],
                            SNP_payload=config["SNPpayload"],
                            Donor_type=config["Donor_type"],
                            Strand_choice=config["Strand_choice"],
                            ssODN_max_size=ssODN_max_size,
                            recoding_args=recoding_args,
                            syn_check_args=syn_check_args,
                            coordinate_without_ENST = coordinate_without_ENST,
                        )
                    except Exception as e:
                        print("Unexpected error:", str(sys.exc_info()))
                        traceback.print_exc()
                        print("additional information:", e)
                        PrintException()
                        guide_seq = current_gRNA["seq"].values[0] if "seq" in current_gRNA.columns else ""
                        csvout_res.write(
                            f"{Entry},{ENST_ID},ERROR: HDR design failed for guide {guide_seq}: {str(e).replace(',', ';')}\n"
                        )
                        continue

                    # append the best gRNA to the final df
                    best_stop_gRNAs = pd.concat([best_stop_gRNAs, current_gRNA])

                    # append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                        cfd1 = HDR_template.cfd_score_post_mut_ins
                    if not hasattr(HDR_template, "cfd_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cfd_score_post_mut2
                    if not hasattr(HDR_template, "cfd_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cfd_score_post_mut3
                    if not hasattr(HDR_template, "cfd_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cfd_score_post_mut4

                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                    # write csv
                    (
                        spec_score,
                        seq,
                        pam,
                        s,
                        e,
                        cut2ins_dist,
                        spec_weight,
                        dist_weight,
                        pos_weight,
                        final_weight,
                    ) = get_res(current_gRNA, spec_score_flavor)

                    donor = HDR_template.Donor_final
                    donor_trimmed = "N/A for ssODN"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed = HDR_template.Donor_final
                        donor = HDR_template.Donor_pretrim

                    # gRNA and donor names
                    # Include gRNA rank in the name when multiple gRNAs are requested
                    gRNA_rank_suffix = f"_rank{i+1}" if gRNA_num_out > 1 else ""

                    ENST_design_counts[ENST_ID] = ENST_design_counts.get(ENST_ID, 0) + 1
                    gRNA_name = f"{ENST_ID}_C_gRNA_entry{Entry}{gRNA_rank_suffix}"
                    donor_name = f"{ENST_ID}_C_donor_entry{Entry}{gRNA_rank_suffix}"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed_name = donor_name.replace("donor_", "donor_trimmed_")
                        #donor_trimmed_name = f"{ENST_ID}_donor_trimmed_{ENST_design_counts[ENST_ID]}"
                    else:
                        donor_trimmed_name = "N/A for ssODN"

                    gRNA_cut_pos = (
                        HDR_template.CutPos
                    )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_C.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        csvout_res.write(
                            f"{Entry},{row_prefix},C,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write( f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )
                    else:
                        csvout_C.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        if not isinstance(cfd4, float):
                            cfd4 = ""
                        csvout_res.write(
                            f"{Entry},{row_prefix},C,{i+1},{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{chopchop_rank},{target_region_label},{target_region_start},{target_region_end},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )

                    # donor features
                    donor_features = HDR_template.Donor_features
                    write_ha_csv_row(ha_out, Entry, ENST_ID, "C", i+1, gRNA_name, seq, insert_pos, HDR_template)
                    # write genbank file
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                        write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])
                    # write anothergenbank file without payload
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}_gRNAonly_noPayload.gb"), "w") as gb_handle:
                        write_genbank_gRNAonly_noPayload(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"], payload_type = config["payload_type"])

                    # write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n   donor before any recoding:{HDR_template.Donor_vanillia}\n    donor after all recoding:{HDR_template.Donor_postMut}\n             donor centered:{HDR_template.Donor_final}\ndonor centered (best strand):{HDR_template.Donor_final}\n\n"
                    #this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template, "info_phase4_5UTR"):
                        fiveUTR_log.write(
                            f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                        )
                    if hasattr(HDR_template, "info_phase5_5UTR"):
                        fiveUTR_log.write(
                            f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                        )

            protein_coding_transcripts_count += 1
            # else:
            #     log.info(f"skipping {ENST_ID} transcript type: {transcript_type} b/c transcript is not protein_coding")
            transcript_count += 1
            # report progress
            if (
                protein_coding_transcripts_count % 100 == 0
                and protein_coding_transcripts_count != 0
            ):
                endtime = datetime.datetime.now()
                elapsed_sec = endtime - starttime
                elapsed_min = elapsed_sec.seconds / 60
                skipped_count = transcript_count - protein_coding_transcripts_count
                log.info(
                    f"processed {protein_coding_transcripts_count}/{total_entries} transcripts (skipped {skipped_count}), elapsed time {elapsed_min:.2f} min ({elapsed_sec} sec)"
                )

        # write csv out
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        if "num_to_process" in locals():
            pass
        else:
            num_to_process = "all"
        skipped_count = max(0, total_entries - protein_coding_transcripts_count)

        log.info(
            f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{total_entries} transcripts\n{skipped_count} nonprotein-coding transcripts were skipped\n"
            f"results written to {config['outdir']}"
        )

        recut_CFD_all.close()
        recut_CFD_fail.close()
        csvout_N.close()
        csvout_C.close()
        fiveUTR_log.close()
        csvout_res.close()
        csvout_res2.close()
        guides_out.close()
        ha_out.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

## end of main()

def translate_sequence(dna_sequence):
    return str(Seq(dna_sequence).translate())


def _valid_coord_pair(v):
    try:
        if v is None:
            return False
        if not isinstance(v, (list, tuple)):
            return False
        if len(v) < 2:
            return False
        int(v[0])
        int(v[1])
        return True
    except Exception:
        return False


def _genbank_record_names(donor_name):
    record_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(donor_name)).strip("._")
    if record_id == "":
        record_id = "protoSpaceJAM_donor"

    locus = record_id
    if len(locus) > 16:
        locus = locus.replace("ENST", "")
        locus = locus.lstrip("0")
        locus = locus.replace("donor", "")
        locus = locus.replace("donor_trimmed", "")
        locus = locus.replace("entry", "")
        locus = locus.replace("__", "_")
        if len(locus) > 16:
            locus = locus[:16]
    if locus == "":
        locus = record_id[:16] if len(record_id) > 16 else record_id
    return locus, record_id

def write_genbank(handle, data_obj, donor_name, donor_type, payload_type):
    """write genebank file"""
    donor_locus, donor_record_id = _genbank_record_names(donor_name)

    gb_record = Record.Record()
    gb_record.locus = donor_locus
    gb_record.size = len(data_obj.Donor_final)
    gb_record.residue_type = 'DNA'
    gb_record.data_file_division = 'PLN'
    gb_record.definition = ''
    gb_record.accession = [donor_record_id]
    gb_record.version = ''
    gb_record.keywords = ['DNA donor']
    gb_record.source = "synthetic DNA donor"
    gb_record.organism = 'synthetic DNA donor'
    gb_record.taxonomy = ['synthetic DNA donor']
    gb_record.sequence=data_obj.Donor_final

    # Convert the Bio.GenBank record to a SeqRecord (needed for writing with SeqIO)
    sequence = Seq(gb_record.sequence)
    seq_record = SeqRecord(sequence, id=donor_record_id, name=donor_locus, description= f"DNA donor type: {donor_type}")
    seq_record.annotations["date"] = get_current_date_formatted()
    seq_record.annotations["data_file_division"] = gb_record.data_file_division
    seq_record.annotations["organism"] = 'synthetic DNA donor'
    seq_record.annotations["molecule_type"] = "DNA"
    seq_record.annotations["accessions"] = [donor_record_id]

    donor_features = getattr(data_obj, "Donor_features", None)
    if not isinstance(donor_features, dict):
        donor_features = {}

    # Features
    if _valid_coord_pair(donor_features.get("left_arm_coord")) and ("HA_payload_strand" in donor_features):
        feature = SeqFeature(FeatureLocation(start=donor_features["left_arm_coord"][0], end=donor_features["left_arm_coord"][1], strand=donor_features["HA_payload_strand"]), type='misc_feature', qualifiers={"label": "left HA"})
        seq_record.features.append(feature)

    if _valid_coord_pair(donor_features.get("right_arm_coord")) and ("HA_payload_strand" in donor_features):
        feature = SeqFeature(FeatureLocation(start=donor_features["right_arm_coord"][0], end=donor_features["right_arm_coord"][1], strand=donor_features["HA_payload_strand"]), type='misc_feature', qualifiers={"label": "right HA"})
        seq_record.features.append(feature)

    if _valid_coord_pair(donor_features.get("tag_coord")) and ("HA_payload_strand" in donor_features):
        feature = SeqFeature(FeatureLocation(start=donor_features["tag_coord"][0], end=donor_features["tag_coord"][1], strand=donor_features["HA_payload_strand"]), type='misc_feature', qualifiers={"label": "payload"})
        seq_record.features.append(feature)

    if "coding_coord" in donor_features and "HA_payload_strand" in donor_features:
        for feat in donor_features["coding_coord"]:
            if not _valid_coord_pair(feat):
                continue
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), strand=donor_features["HA_payload_strand"], type='exon', qualifiers={"label": "exon"})
            seq_record.features.append(feature)

    if "gRNA_coord" in donor_features and "gRNA_strand" in donor_features:
        for feat in donor_features["gRNA_coord"]:
            if not _valid_coord_pair(feat):
                continue
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1], strand=donor_features["gRNA_strand"]), type='misc_feature', qualifiers={"label": "gRNA+PAM"})
            seq_record.features.append(feature)
    if "recoding_coord" in donor_features:
        for feat in donor_features["recoding_coord"]:
            if not _valid_coord_pair(feat):
                continue
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), type='misc_feature', qualifiers={"label": "recode"})
            seq_record.features.append(feature)

    # Write to a GenBank file using SeqIO for the actual file writing
    SeqIO.write([seq_record], handle, 'genbank')

def write_genbank_gRNAonly_noPayload(handle, data_obj, donor_name, donor_type, payload_type):
    """write genebank file with wt homology arm and gRNA only"""
    donor_locus, donor_record_id = _genbank_record_names(donor_name)

    gb_record = Record.Record()
    gb_record.locus = donor_locus
    gb_record.size = len(data_obj.left_flk_seq + data_obj.right_flk_seq)
    gb_record.residue_type = 'DNA'
    gb_record.data_file_division = 'PLN'
    gb_record.definition = ''
    gb_record.accession = [donor_record_id]
    gb_record.version = ''
    gb_record.keywords = ['DNA sequence']
    gb_record.source = "DNA sequence"
    gb_record.organism = 'DNA sequence'
    gb_record.taxonomy = ['DNA sequence']
    gb_record.sequence=data_obj.left_flk_seq + data_obj.right_flk_seq

    # Convert the Bio.GenBank record to a SeqRecord (needed for writing with SeqIO)
    sequence = Seq(gb_record.sequence)
    seq_record = SeqRecord(sequence, id=donor_record_id, name=donor_locus, description= f"wt homology arms and gRNA")
    seq_record.annotations["date"] = get_current_date_formatted()
    seq_record.annotations["data_file_division"] = gb_record.data_file_division
    seq_record.annotations["organism"] = 'DNA sequence'
    seq_record.annotations["molecule_type"] = "DNA"
    seq_record.annotations["accessions"] = [donor_record_id]

    payloadless_features = getattr(data_obj, "payloadless_donor_features", None)
    if not isinstance(payloadless_features, dict):
        payloadless_features = {}

    # Features
    if _valid_coord_pair(payloadless_features.get("left_arm_coord")) and ("HA_payload_strand" in payloadless_features):
        feature = SeqFeature(FeatureLocation(start=payloadless_features["left_arm_coord"][0], end=payloadless_features["left_arm_coord"][1], strand=payloadless_features["HA_payload_strand"]), type='misc_feature', qualifiers={"label": "left homology arm (before trimming)"})
        seq_record.features.append(feature)

    if _valid_coord_pair(payloadless_features.get("right_arm_coord")) and ("HA_payload_strand" in payloadless_features):
        feature = SeqFeature(FeatureLocation(start=payloadless_features["right_arm_coord"][0], end=payloadless_features["right_arm_coord"][1], strand=payloadless_features["HA_payload_strand"]), type='misc_feature', qualifiers={"label": "right homology arm (before trimming)"})
        seq_record.features.append(feature)

    if "coding_coord" in payloadless_features and "HA_payload_strand" in payloadless_features:
        for feat in payloadless_features["coding_coord"]:
            if not _valid_coord_pair(feat):
                continue
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), strand=payloadless_features["HA_payload_strand"], type='exon', qualifiers={"label": "exon"})
            seq_record.features.append(feature)

    if "ORF_coord" in payloadless_features and "HA_payload_strand" in payloadless_features:
        for feat in payloadless_features["ORF_coord"]:
            if not _valid_coord_pair(feat):
                continue
            #feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), strand=data_obj.Donor_features["HA_payload_strand"], type='CDS-in-frame', qualifiers={"label": "CDS-in-frame"})
            #seq_record.features.append(feature)

            # Extract and translate the coding in-frame sequence
            orf_sequence = sequence[feat[0]:feat[1]]
            protein_sequence = translate_sequence(str(orf_sequence))
            protein_feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), strand=payloadless_features["HA_payload_strand"],type='CDS', qualifiers={"label": "CDS", "codon_start": 1, "translation": protein_sequence})
            seq_record.features.append(protein_feature)

    if "gRNA_coord" in payloadless_features and "gRNA_strand" in payloadless_features:
        for feat in payloadless_features["gRNA_coord"]:
            if not _valid_coord_pair(feat):
                continue
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1], strand=payloadless_features["gRNA_strand"]), type='gRNA+PAM', qualifiers={"label": "gRNA+PAM"})
            seq_record.features.append(feature)


    # Write to a GenBank file using SeqIO for the actual file writing
    SeqIO.write([seq_record], handle, 'genbank')


def ret_six_dec(myvar):
    """
    retain six decimal points for printout
    """
    if type(myvar) == float:
        return f"{myvar:.6f}"
    elif type(myvar) == int:
        myvar = float(myvar)
        return f"{myvar:.6f}"
    else:
        return myvar


def mkdir(mypath):
    if not os.path.exists(mypath):
        os.makedirs(mypath)


def deepmerge(dict1, dict2):
    """
    merge two dictionary at the secondary key level
    """
    if len(dict1) == 0:
        return dict2
    if len(dict2) == 0:
        return dict1
    # start merging
    dictm = dict1
    for k in dict2:
        if k in dict1:  # shared key
            for k2 in dict2[k].keys():
                dictm[k][k2] = dict2[k][k2]
        else:
            dictm[k] = dict2[k]
    return dictm


class info:
    """
    info log class
    """
    def __init__(self) -> None:
        self.cfd1 = []
        self.cfd2 = []
        self.cfd3 = []
        self.cfd4 = []
        self.cfdfinal = []
        self.failed = []


def get_res(best_start_gRNA, spec_score_flavor):
    spec_score = best_start_gRNA[spec_score_flavor].values[0]
    seq = best_start_gRNA["seq"].values[0]
    pam = best_start_gRNA["pam"].values[0]
    s = best_start_gRNA["start"].values[0]
    e = best_start_gRNA["end"].values[0]
    cut2ins_dist = best_start_gRNA["Cut2Ins_dist"].values[0]
    spec_weight = best_start_gRNA["spec_weight"].values[0]
    dist_weight = best_start_gRNA["dist_weight"].values[0]
    pos_weight = best_start_gRNA["pos_weight"].values[0]
    final_weight = best_start_gRNA["final_weight"].values[0]
    return [
        spec_score,
        seq,
        pam,
        s,
        e,
        cut2ins_dist,
        spec_weight,
        dist_weight,
        pos_weight,
        final_weight,
    ]


def load_sequence_from_file(path):
    """
    Load a nucleotide sequence from plain text or FASTA.
    Removes whitespace and ignores FASTA header lines.
    """
    if not os.path.isfile(path):
        sys.exit(f"payload file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    seq_parts = []
    for raw in lines:
        line = raw.strip()
        if line == "":
            continue
        if line.startswith(">"):
            continue
        seq_parts.append(line)
    seq = "".join(seq_parts).replace(" ", "").replace("\t", "")
    if seq == "":
        sys.exit(f"payload file is empty or contains no sequence: {path}")
    return seq


def resolve_payload_inputs_from_files(config):
    """
    Resolve optional payload file inputs and write back into config payload fields.
    File-based values override direct string payload args.
    """
    mapping = [
        ("payload_file", "payload"),
        ("Npayload_file", "Npayload"),
        ("Cpayload_file", "Cpayload"),
        ("POSpayload_file", "POSpayload"),
        ("SNPpayload_file", "SNPpayload"),
    ]
    for file_key, seq_key in mapping:
        fpath = str(config.get(file_key, "")).strip()
        if fpath != "":
            config[seq_key] = load_sequence_from_file(fpath)


def validate_payload_inputs(config, log):
    payload_type = str(config.get("payload_type", "insertion")).strip().lower()
    if payload_type != "insertion":
        return

    payload_specs = [
        ("payload", "global payload"),
        ("Npayload", "N-terminus payload"),
        ("Cpayload", "C-terminus payload"),
        ("POSpayload", "coordinate payload"),
    ]
    for key, label in payload_specs:
        seq = str(config.get(key, "") or "").strip().upper()
        if seq == "":
            continue
        clean = re.sub(r"[^ACGTN]", "", seq)
        if clean == "":
            log.warning(f"{label} is empty after removing non-ACGTN characters")
            continue
        if clean != seq:
            log.warning(
                f"{label} contains non-ACGTN characters; frame check uses sanitized length {len(clean)} bp"
            )
        mod = len(clean) % 3
        if mod == 0:
            log.info(f"{label}: {len(clean)} bp (in-frame, length mod 3 = 0)")
        else:
            log.warning(f"{label}: {len(clean)} bp (out-of-frame, length mod 3 = {mod})")


def load_input_dataframe(config, outdir, log):
    """Load target table from CSV or build one-row table from direct CLI args."""
    if str(config["input_mode"]).lower() == "direct":
        direct_ensembl_id = str(config["direct_ensembl_id"]).strip()
        direct_chr = str(config["direct_chromosome"]).strip()
        direct_coord = str(config["direct_coordinate"]).strip()
        direct_term = str(config["direct_target_terminus"]).strip().upper()
        direct_gene_name = str(config.get("direct_gene_name", "")).strip()

        if direct_term not in ["N", "C", "ALL", ""]:
            sys.exit("--direct_target_terminus must be N, C, ALL, or empty")
        if direct_coord != "" and not direct_coord.isdigit():
            sys.exit("--direct_coordinate must be an integer")

        has_coordinate = (direct_chr != "" and direct_coord != "")
        has_enst = (direct_ensembl_id != "")
        has_gene = (direct_gene_name != "")
        if (not has_coordinate) and (not has_enst) and (not has_gene):
            sys.exit(
                "For --input_mode direct, provide either --direct_ensembl_id (ENST mode) "
                "or both --direct_chromosome and --direct_coordinate (coordinate mode) "
                "or --direct_gene_name (CHOPCHOP web gene query mode)"
            )

        row = {
            "Entry": str(config["direct_entry"]),
            "Ensembl_ID": direct_ensembl_id,
            "Target_terminus": direct_term if direct_term != "" else "ALL",
        }
        if has_coordinate:
            row["Chromosome"] = direct_chr
            row["Coordinate"] = direct_coord
        if has_gene:
            row["Gene_Name"] = direct_gene_name
        df = pd.DataFrame([row], dtype=str)
        auto_input_csv = os.path.join(outdir, "auto_input_from_cli.csv")
        df.to_csv(auto_input_csv, index=False)
        log.info(
            f"input_mode=direct: created one-row input table from CLI args ({auto_input_csv})"
        )
        return df

    if os.path.isfile(config["path2csv"]):
        log.info(
            f"begin processing user-supplied list of gene IDs in file {config['path2csv']}"
        )
        df = pd.read_csv(os.path.join(config["path2csv"]), dtype=str)
        keys2check = set(["Ensembl_ID"])
        if not keys2check.issubset(df.columns):
            log.error(
                'Missing columns in the input csv file\n Required columns:"Ensembl_ID"'
            )
            log.info("Please fix the input csv file and try again")
            sys.exit()
        return df

    sys.exit(f"ERROR: The input file {config['path2csv']} is not found")


def load_guides_dataframe(config, log):
    path = str(config.get("guides_csv", "")).strip()
    if path == "":
        return None
    if not os.path.isfile(path):
        sys.exit(f"ERROR: guides_csv file not found: {path}")

    df = pd.read_csv(path, dtype=str)
    rename_map = {
        "insert_pos": "Insert_pos",
        "cut2ins_dist": "Cut2Ins_dist",
        "eff_scores": "Eff_scores",
    }
    df = df.rename(columns=rename_map)
    for col in ["start", "end", "Insert_pos", "Cut2Ins_dist"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Ensure required columns exist for downstream HDR and reporting
    defaults = {
        "seq": "",
        "pam": "",
        "start": 0,
        "end": 0,
        "strand": "+",
        "guideMITScore": 50,
        "guideCfdScore": 50,
        "guideCfdScorev2": 50,
        "guideCfdScorev3": 50,
        "Eff_scores": 0,
        "chopchop_offtarget_penalty": float("nan"),
        "chr": "",
        "ID": "",
        "Insert_pos": 0,
        "Cut2Ins_dist": 0,
        "spec_weight": 1,
        "dist_weight": 1,
        "pos_weight": 1,
        "final_weight": 1,
        "score_note": "",
    }
    for col, value in defaults.items():
        if col not in df.columns:
            df[col] = value
        df[col] = df[col].fillna(value)

    if "rank" in df.columns:
        df["_rank_sort"] = pd.to_numeric(df["rank"], errors="coerce")
        df = df.sort_values("_rank_sort", kind="stable", na_position="last")
    df = df.reset_index(drop=True)
    log.info(f"loaded guides_csv with {df.shape[0]} guides from {path}")
    return df


GENOME_TO_ENSEMBL_SPECIES = {
    "GRCh38": "homo_sapiens",
    "GRCm39": "mus_musculus",
    "GRCz11": "danio_rerio",
    "mRatBN7.2": "rattus_norvegicus",
}


def _ensembl_species(genome_ver):
    return GENOME_TO_ENSEMBL_SPECIES.get(str(genome_ver), "homo_sapiens")


def _ensembl_json(url, timeout=30):
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as exc:
        # Fallback for environments with broken cert bundle.
        if "CERTIFICATE_VERIFY_FAILED" in str(exc) or isinstance(exc, ssl.SSLError):
            ctx = ssl._create_unverified_context()
            with urlopen(req, timeout=timeout, context=ctx) as response:
                return json.loads(response.read().decode("utf-8"))
        raise


def _ensembl_json_any(urls, timeout=30):
    last_exc = None
    for url in urls:
        try:
            return _ensembl_json(url, timeout=timeout)
        except Exception as exc:
            last_exc = RuntimeError(f"{exc} @ {url}")
    if last_exc is not None:
        raise last_exc
    raise RuntimeError("No Ensembl URL candidates provided")


def _normalize_interval(start, end):
    a = int(start)
    b = int(end)
    return (a, b) if a <= b else (b, a)


def _normalize_enst_id(enst_id):
    raw = str(enst_id).strip()
    m = re.search(r"(ENST\d+)", raw)
    if m is not None:
        raw = m.group(1)
    if "." in raw:
        raw = raw.split(".")[0]
    return raw


def _add_loc_type(loc_map, start, end, label):
    s, e = _normalize_interval(start, end)
    loc_map[(s, e)] = label


def _get_cache_path(cache_dir, genome_ver, enst_id):
    if not cache_dir:
        return ""
    safe_enst = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(enst_id))
    safe_genome = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(genome_ver))
    return os.path.join(cache_dir, f"{safe_genome}_{safe_enst}.annotation.json")


def _get_regions_cache_path(cache_dir, genome_ver, enst_id):
    if not cache_dir:
        return ""
    safe_enst = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(enst_id))
    safe_genome = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(genome_ver))
    return os.path.join(cache_dir, f"{safe_genome}_{safe_enst}.regions.json")


def _build_transcript_obj_from_bundle_data(data):
    enst_id = data["enst_id"]
    chr_name = str(data["chr"])
    strand = int(data["strand"])
    features = [SimpleNamespace(strand=strand, type="transcript")]

    cds_for_features = data.get("cds", [])
    if strand == 1:
        cds_for_features = sorted(cds_for_features, key=lambda x: int(x[0]))
    else:
        cds_for_features = sorted(cds_for_features, key=lambda x: int(x[0]), reverse=True)
    for s, e in cds_for_features:
        features.append(
            SimpleNamespace(
                type="CDS",
                strand=strand,
                location=SimpleNamespace(ref=chr_name, start=int(s), end=int(e)),
            )
        )
    for s, e in data.get("exons", []):
        features.append(
            SimpleNamespace(
                type="exon",
                strand=strand,
                location=SimpleNamespace(ref=chr_name, start=int(s), end=int(e)),
            )
        )

    transcript_obj = SimpleNamespace(
        chr=chr_name,
        span_start=int(data["span_start"]),
        span_end=int(data["span_end"]),
        name=str(data.get("gene_name", enst_id)),
        description=f"{enst_id}|{str(data.get('biotype', 'unknown'))}|ensembl_rest",
        features=features,
    )
    return {enst_id: transcript_obj}


def _region_record(label, intervals, strand, region_type="", exon_number=None):
    normalized = [_normalize_interval(s, e) for s, e in intervals]
    if len(normalized) == 0:
        return None
    positions = _positions_from_intervals(normalized, strand)
    return {
        "label": str(label),
        "region_type": str(region_type) if region_type else str(label),
        "exon_number": int(exon_number) if exon_number is not None else None,
        "intervals": [[int(s), int(e)] for s, e in normalized],
        "start": int(min(positions)),
        "end": int(max(positions)),
        "transcript_anchor_start": int(positions[0]),
        "transcript_anchor_end": int(positions[-1]),
    }


def _build_region_catalog(chr_name, strand, exons, cds, utr5, utr3):
    strand = int(strand)
    exon_order = _intervals_in_transcript_order(exons, strand)
    regions = []

    def add_region(label, intervals, region_type="", exon_number=None):
        record = _region_record(
            label=label,
            intervals=intervals,
            strand=strand,
            region_type=region_type,
            exon_number=exon_number,
        )
        if record is not None:
            regions.append(record)

    add_region("transcript", exons, "transcript")
    add_region("CDS", cds, "CDS")
    add_region("5UTR", utr5, "5UTR")
    add_region("3UTR", utr3, "3UTR")

    for exon_num, exon_interval in enumerate(exon_order, start=1):
        add_region(f"exon{exon_num}", [exon_interval], "exon", exon_num)
        overlaps = []
        for cds_start, cds_end in cds:
            ov_start = max(int(exon_interval[0]), int(cds_start))
            ov_end = min(int(exon_interval[1]), int(cds_end))
            if ov_start <= ov_end:
                overlaps.append([ov_start, ov_end])
        if len(overlaps) > 0:
            add_region(f"CDS_exon{exon_num}", overlaps, "CDS", exon_num)

    return {
        "chr": str(chr_name),
        "strand": strand,
        "regions": regions,
    }


def _build_region_catalog_from_bundle_data(bundle_data):
    return _build_region_catalog(
        chr_name=bundle_data.get("chr", ""),
        strand=bundle_data.get("strand", 1),
        exons=bundle_data.get("exons", []),
        cds=bundle_data.get("cds", []),
        utr5=bundle_data.get("utr5", []),
        utr3=bundle_data.get("utr3", []),
    )


def _build_region_catalog_from_transcript_obj(transcript_obj):
    strand = int(getattr(transcript_obj.features[0], "strand", 1))
    exons = _transcript_feature_intervals(transcript_obj, "exon")
    cds = _transcript_feature_intervals(transcript_obj, "CDS")
    utr5, utr3 = _derive_utr_from_exons_and_cds(exons=exons, cds=cds, strand=strand)
    return _build_region_catalog(
        chr_name=str(getattr(transcript_obj, "chr", "")),
        strand=strand,
        exons=exons,
        cds=cds,
        utr5=utr5,
        utr3=utr3,
    )


def _get_region_entry_from_catalog(region_catalog, region="", exon_num=None):
    requested = _normalize_preferred_region(region)
    if requested == "":
        requested = "CDS"
    if requested == "exon":
        if exon_num is None or exon_num <= 0:
            raise ValueError("Preferred_exon is required and must be >= 1 when Preferred_region=exon")
        target_label = f"exon{int(exon_num)}"
    elif requested == "CDS":
        target_label = f"CDS_exon{int(exon_num)}" if exon_num is not None else "CDS"
    elif requested in ["5UTR", "3UTR", "transcript"]:
        target_label = requested
    else:
        target_label = "transcript"

    for region_entry in region_catalog.get("regions", []):
        if str(region_entry.get("label", "")) == target_label:
            return region_entry
    if requested == "CDS" and exon_num is not None:
        raise ValueError("Preferred_exon=%s has no CDS segment in this transcript" % exon_num)
    raise ValueError("Requested region '%s' is unavailable for this transcript" % target_label)


def _build_loc2postype_for_transcript(enst_id, chr_name, strand, exons, cds, utr5, utr3):
    mapping = {}
    for s, e in cds:
        _add_loc_type(mapping, s, e, "cds")
    for s, e in utr5:
        _add_loc_type(mapping, s, e, "5UTR")
    for s, e in utr3:
        _add_loc_type(mapping, s, e, "3UTR")

    exons_sorted = sorted([_normalize_interval(s, e) for s, e in exons], key=lambda x: x[0])
    if len(exons_sorted) >= 2:
        if int(strand) == 1:
            exon_order = exons_sorted
        else:
            exon_order = list(reversed(exons_sorted))

        for i in range(len(exon_order) - 1):
            curr_exon = exon_order[i]
            next_exon = exon_order[i + 1]
            curr_start, curr_end = curr_exon
            next_start, next_end = next_exon

            if int(strand) == 1:
                exin_j = curr_end
                inex_j = next_start
                _add_loc_type(mapping, exin_j - 1, exin_j + 2, "within_2bp_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 2, exin_j + 3, "within_3bp_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 3, exin_j - 2, "3N4bp_up_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j + 3, exin_j + 6, "3_to_6bp_down_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 2, exin_j + 6, "-3_to_+6bp_of_exon_intron_junction")

                _add_loc_type(mapping, inex_j - 2, inex_j + 1, "within_2bp_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 3, inex_j + 2, "within_3bp_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 4, inex_j - 3, "3N4bp_up_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j + 2, inex_j + 3, "3N4bp_down_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 3, inex_j + 1, "-3_to_+2bp_of_intron_exon_junction")
            else:
                exin_j = curr_start
                inex_j = next_end
                _add_loc_type(mapping, exin_j - 2, exin_j + 1, "within_2bp_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 3, exin_j + 2, "within_3bp_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j + 2, exin_j + 3, "3N4bp_up_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 6, exin_j - 3, "3_to_6bp_down_of_exon_intron_junction")
                _add_loc_type(mapping, exin_j - 6, exin_j + 2, "-3_to_+6bp_of_exon_intron_junction")

                _add_loc_type(mapping, inex_j - 1, inex_j + 2, "within_2bp_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 2, inex_j + 3, "within_3bp_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j + 3, inex_j + 4, "3N4bp_up_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 3, inex_j - 2, "3N4bp_down_of_intron_exon_junction")
                _add_loc_type(mapping, inex_j - 1, inex_j + 3, "-3_to_+2bp_of_intron_exon_junction")

    return {str(chr_name): {str(enst_id): mapping}}


def _build_phase_map_for_transcript(enst_id, chr_name, strand, cds):
    phase_map = {}
    cds_norm = [_normalize_interval(s, e) for s, e in cds]
    if len(cds_norm) == 0:
        return {}
    if int(strand) == 1:
        cds_order = sorted(cds_norm, key=lambda x: x[0])
    else:
        cds_order = sorted(cds_norm, key=lambda x: x[0], reverse=True)

    coding_positions = []
    for s, e in cds_order:
        if int(strand) == 1:
            coding_positions.extend(range(int(s), int(e) + 1))
        else:
            coding_positions.extend(range(int(e), int(s) - 1, -1))

    chr_map = {}
    for idx, pos in enumerate(coding_positions):
        phase = (idx % 3) + 1
        if pos not in chr_map:
            chr_map[pos] = {}
        chr_map[pos][str(enst_id)] = phase
    if len(chr_map) == 0:
        return {}
    phase_map[str(chr_name)] = chr_map
    return phase_map


def _derive_utr_from_exons_and_cds(exons, cds, strand):
    utr5 = []
    utr3 = []
    exons_norm = [_normalize_interval(s, e) for s, e in exons]
    cds_norm = [_normalize_interval(s, e) for s, e in cds]
    if len(exons_norm) == 0 or len(cds_norm) == 0:
        return utr5, utr3

    cds_min = min([s for s, _ in cds_norm])
    cds_max = max([e for _, e in cds_norm])
    for ex_s, ex_e in exons_norm:
        if int(strand) == 1:
            if ex_s < cds_min:
                left_s = ex_s
                left_e = min(ex_e, cds_min - 1)
                if left_s <= left_e:
                    utr5.append([left_s, left_e])
            if ex_e > cds_max:
                right_s = max(ex_s, cds_max + 1)
                right_e = ex_e
                if right_s <= right_e:
                    utr3.append([right_s, right_e])
        else:
            if ex_e > cds_max:
                right_s = max(ex_s, cds_max + 1)
                right_e = ex_e
                if right_s <= right_e:
                    utr5.append([right_s, right_e])
            if ex_s < cds_min:
                left_s = ex_s
                left_e = min(ex_e, cds_min - 1)
                if left_s <= left_e:
                    utr3.append([left_s, left_e])
    return utr5, utr3


def get_ensembl_annotation_bundle(enst_id, genome_ver, cache=None, cache_dir="", write_cache=False, timeout=30):
    if not enst_id:
        logging.getLogger("protoSpaceJAM").warning("Ensembl annotation lookup skipped: empty ENST ID")
        return None
    enst_id = _normalize_enst_id(enst_id)
    if cache is not None and enst_id in cache:
        return cache[enst_id]

    bundle_data = None
    cache_path = _get_cache_path(cache_dir, genome_ver, enst_id)
    if cache_path and os.path.isfile(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                bundle_data = json.load(fh)
        except Exception:
            logging.getLogger("protoSpaceJAM").warning(
                f"failed reading Ensembl annotation cache file, refetching: {cache_path}"
            )
            bundle_data = None

    if bundle_data is None:
        try:
            ensembl_bases = ["https://rest.ensembl.org"]
            lookup_urls = [f"{base}/lookup/id/{enst_id}?expand=1" for base in ensembl_bases]
            lookup = _ensembl_json_any(lookup_urls, timeout=timeout)
            if not isinstance(lookup, dict) or "seq_region_name" not in lookup:
                logging.getLogger("protoSpaceJAM").warning(
                    f"Ensembl lookup returned unexpected payload for {enst_id}"
                )
                return None
            chr_name = str(lookup.get("seq_region_name", ""))
            start = int(lookup.get("start", 0))
            end = int(lookup.get("end", 0))
            strand = int(lookup.get("strand", 1))
            biotype = str(lookup.get("biotype", "unknown"))
            gene_name = str(lookup.get("display_name") or lookup.get("external_name") or enst_id)

            overlap_urls = [f"{base}/overlap/id/{enst_id}?feature=exon&feature=cds" for base in ensembl_bases]
            overlap = _ensembl_json_any(overlap_urls, timeout=timeout)
            if not isinstance(overlap, list):
                overlap = []
            exons = []
            cds = []
            for feat in overlap:
                if not isinstance(feat, dict):
                    continue
                fs = int(feat.get("start", 0))
                fe = int(feat.get("end", 0))
                ftype = str(feat.get("feature_type", "")).lower()
                if fs <= 0 or fe <= 0:
                    continue
                if ftype == "exon":
                    exons.append([fs, fe])
                elif ftype == "cds":
                    cds.append([fs, fe])
            utr5, utr3 = _derive_utr_from_exons_and_cds(exons=exons, cds=cds, strand=strand)

            bundle_data = {
                "enst_id": str(enst_id),
                "chr": chr_name,
                "span_start": start,
                "span_end": end,
                "strand": strand,
                "biotype": biotype,
                "gene_name": gene_name,
                "exons": exons,
                "cds": cds,
                "utr5": utr5,
                "utr3": utr3,
            }
        except (HTTPError, URLError, TimeoutError, ValueError, ssl.SSLError, RuntimeError) as exc:
            logging.getLogger("protoSpaceJAM").warning(
                f"Ensembl annotation fetch failed for {enst_id}: {exc}. Check internet/SSL access to Ensembl REST."
            )
            return None

        if write_cache and cache_path:
            try:
                mkdir(os.path.dirname(cache_path))
                with open(cache_path, "w", encoding="utf-8") as fh:
                    json.dump(bundle_data, fh, indent=2)
            except Exception:
                pass

    region_catalog = _build_region_catalog_from_bundle_data(bundle_data)
    regions_cache_path = _get_regions_cache_path(cache_dir, genome_ver, enst_id)
    if regions_cache_path:
        try:
            mkdir(os.path.dirname(regions_cache_path))
            region_payload = {
                "enst_id": str(bundle_data.get("enst_id", enst_id)),
                "gene_name": str(bundle_data.get("gene_name", enst_id)),
                "chr": str(bundle_data.get("chr", "")),
                "strand": int(bundle_data.get("strand", 1)),
                "span_start": int(bundle_data.get("span_start", 0)),
                "span_end": int(bundle_data.get("span_end", 0)),
                "regions": region_catalog.get("regions", []),
            }
            with open(regions_cache_path, "w", encoding="utf-8") as fh:
                json.dump(region_payload, fh, indent=2)
        except Exception:
            pass

    transcript_dict = _build_transcript_obj_from_bundle_data(bundle_data)
    enst_key = str(bundle_data["enst_id"])
    chr_name = str(bundle_data["chr"])
    strand = int(bundle_data["strand"])
    cds = bundle_data.get("cds", [])
    exons = bundle_data.get("exons", [])
    utr5 = bundle_data.get("utr5", [])
    utr3 = bundle_data.get("utr3", [])
    loc2posType = _build_loc2postype_for_transcript(
        enst_id=enst_key,
        chr_name=chr_name,
        strand=strand,
        exons=exons,
        cds=cds,
        utr5=utr5,
        utr3=utr3,
    )
    phase_map = _build_phase_map_for_transcript(
        enst_id=enst_key,
        chr_name=chr_name,
        strand=strand,
        cds=cds,
    )

    out = {
        "ENST_info": transcript_dict,
        "loc2posType": loc2posType,
        "ENST_PhaseInCodon": phase_map,
        "region_catalog": region_catalog,
    }
    if cache is not None:
        cache[enst_id] = out
    return out


def get_insert_positions_from_ensembl_bundle(genome_ver, enst_id, cache=None, cache_dir="", write_cache=False):
    enst_norm = _normalize_enst_id(enst_id)
    bundle = get_ensembl_annotation_bundle(
        enst_id=enst_norm,
        genome_ver=genome_ver,
        cache=cache,
        cache_dir=cache_dir,
        write_cache=write_cache,
    )
    if bundle is None or "ENST_info" not in bundle or enst_norm not in bundle["ENST_info"]:
        return None
    try:
        ENST_info = bundle["ENST_info"]
        ATG_loc, stop_loc = get_start_stop_loc(enst_norm, ENST_info)
        n_insert = int(get_end_pos_of_ATG(ATG_loc)[1])
        c_insert = int(get_start_pos_of_stop(stop_loc)[1])
        chr_name = str(ENST_info[enst_norm].chr)
        return (n_insert, c_insert, chr_name)
    except Exception:
        return None


def _clean_optional_str(value):
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _normalize_preferred_region(value):
    raw = _clean_optional_str(value)
    if raw == "":
        return ""
    compact = re.sub(r"[^a-z0-9]+", "", raw.lower())
    region_map = {
        "exon": "exon",
        "5utr": "5UTR",
        "utr5": "5UTR",
        "fiveutr": "5UTR",
        "fiveprimeutr": "5UTR",
        "5primeutr": "5UTR",
        "3utr": "3UTR",
        "utr3": "3UTR",
        "threeutr": "3UTR",
        "threeprimeutr": "3UTR",
        "3primeutr": "3UTR",
        "cds": "CDS",
        "coding": "CDS",
        "codingsequence": "CDS",
        "transcript": "transcript",
        "mrna": "transcript",
    }
    return region_map.get(compact, "")


def _normalize_preferred_anchor(value):
    raw = _clean_optional_str(value)
    if raw == "":
        return "center"
    compact = re.sub(r"[^a-z0-9]+", "", raw.lower())
    anchor_map = {
        "start": "start",
        "5prime": "start",
        "5p": "start",
        "beginning": "start",
        "end": "end",
        "3prime": "end",
        "3p": "end",
        "stop": "end",
        "center": "center",
        "centre": "center",
        "middle": "center",
        "mid": "center",
    }
    return anchor_map.get(compact, "center")


def _parse_optional_int(value, field_name):
    text = _clean_optional_str(value)
    if text == "":
        return None
    try:
        return int(text)
    except Exception:
        raise ValueError("%s must be an integer if provided" % field_name)


def _intervals_in_transcript_order(intervals, strand):
    normalized = [_normalize_interval(s, e) for s, e in intervals]
    if int(strand) == 1:
        return sorted(normalized, key=lambda x: x[0])
    return sorted(normalized, key=lambda x: x[0], reverse=True)


def _positions_from_intervals(intervals, strand):
    ordered = _intervals_in_transcript_order(intervals, strand)
    positions = []
    for start, end in ordered:
        if int(strand) == 1:
            positions.extend(range(int(start), int(end) + 1))
        else:
            positions.extend(range(int(end), int(start) - 1, -1))
    return positions


def _transcript_feature_intervals(transcript_obj, feature_type):
    wanted = str(feature_type).lower()
    intervals = []
    for feat in getattr(transcript_obj, "features", []):
        feat_type = str(getattr(feat, "type", "")).lower()
        if feat_type != wanted:
            continue
        location = getattr(feat, "location", None)
        if location is None:
            continue
        start = getattr(location, "start", None)
        end = getattr(location, "end", None)
        if start is None or end is None:
            continue
        intervals.append([int(start), int(end)])
    return intervals


def _resolve_anchor_coordinate_from_positions(positions, anchor="center", offset=0):
    if len(positions) == 0:
        raise ValueError("No genomic positions available for the requested region")
    anchor = _normalize_preferred_anchor(anchor)
    if anchor == "start":
        idx = 0
    elif anchor == "end":
        idx = len(positions) - 1
    else:
        idx = (len(positions) - 1) // 2
    idx = idx + int(offset)
    idx = max(0, min(len(positions) - 1, idx))
    return int(positions[idx])


def resolve_preferred_insert_site(ENST_info, ENST_ID, preferred_region="", preferred_exon="", preferred_anchor="", preferred_offset=""):
    region = _normalize_preferred_region(preferred_region)
    if region == "":
        region = "CDS"
    if ENST_ID not in ENST_info:
        raise ValueError("Transcript annotation is unavailable for %s" % ENST_ID)

    transcript_obj = ENST_info[ENST_ID]
    strand = int(getattr(transcript_obj.features[0], "strand", 1))
    chr_name = str(getattr(transcript_obj, "chr", ""))
    anchor = _normalize_preferred_anchor(preferred_anchor)
    offset = _parse_optional_int(preferred_offset, "Preferred_offset") or 0
    exon_num = _parse_optional_int(preferred_exon, "Preferred_exon")

    region_catalog = _build_region_catalog_from_transcript_obj(transcript_obj)
    region_entry = _get_region_entry_from_catalog(
        region_catalog=region_catalog,
        region=region,
        exon_num=exon_num,
    )

    positions = _positions_from_intervals(region_entry.get("intervals", []), strand)
    coordinate = _resolve_anchor_coordinate_from_positions(
        positions=positions,
        anchor=anchor,
        offset=offset,
    )
    return {
        "chrom": chr_name,
        "coordinate": coordinate,
        "region_label": str(region_entry.get("label", region)),
        "region_start": int(region_entry.get("start")),
        "region_end": int(region_entry.get("end")),
        "anchor": anchor,
        "offset": offset,
    }


def _extract_donor_feature_seq(donor_seq, coord_pair):
    if donor_seq is None or not _valid_coord_pair(coord_pair):
        return ""
    seq_len = len(donor_seq)
    start = max(0, int(coord_pair[0]))
    end = min(seq_len, int(coord_pair[1]))
    if end <= start and end < seq_len:
        end = min(seq_len, end + 1)
    if end <= start:
        return ""
    return str(donor_seq[start:end])


def write_ha_csv_row(handle, entry, enst_id, terminus, design_rank, gRNA_name, gRNA_seq, insert_pos, hdr_template):
    donor_seq = str(getattr(hdr_template, "Donor_final", "") or "")
    donor_features = getattr(hdr_template, "Donor_features", {}) or {}
    left_ha = _extract_donor_feature_seq(donor_seq, donor_features.get("left_arm_coord"))
    right_ha = _extract_donor_feature_seq(donor_seq, donor_features.get("right_arm_coord"))
    payload_seq = _extract_donor_feature_seq(donor_seq, donor_features.get("tag_coord"))
    handle.write(
        "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            entry,
            enst_id,
            terminus,
            design_rank,
            gRNA_name,
            gRNA_seq,
            insert_pos,
            left_ha,
            payload_seq,
            right_ha,
            donor_seq,
        )
    )


def get_insert_positions_from_local_enst(genome_ver, enst_id, cache=None):
    """
    Resolve N- and C-terminus insert positions for an ENST using local parsed annotations.
    Returns (n_insert_pos, c_insert_pos, chr) or None if unavailable.
    """
    if not enst_id:
        return None
    if cache is not None and enst_id in cache:
        return cache[enst_id]

    try:
        idx_path = os.path.join(
            "genome_files",
            "parsed_gff3",
            genome_ver,
            "ENST_info",
            "ENST_info_index.pickle",
        )
        if not os.path.isfile(idx_path):
            return None
        with open(idx_path, "rb") as fh:
            enst_info_index = pickle.load(fh)
        if enst_id not in enst_info_index:
            return None
        part = enst_info_index[enst_id]
        info_path = os.path.join(
            "genome_files",
            "parsed_gff3",
            genome_ver,
            "ENST_info",
            f"ENST_info_part{part}.pickle",
        )
        if not os.path.isfile(info_path):
            return None
        with open(info_path, "rb") as fh:
            ENST_info = pickle.load(fh)
        if enst_id not in ENST_info:
            return None
        ATG_loc, stop_loc = get_start_stop_loc(enst_id, ENST_info)
        n_insert = int(get_end_pos_of_ATG(ATG_loc)[1])
        c_insert = int(get_start_pos_of_stop(stop_loc)[1])
        chr_name = str(ENST_info[enst_id].chr)
        out = (n_insert, c_insert, chr_name)
        if cache is not None:
            cache[enst_id] = out
        return out
    except Exception:
        return None


def _cut_pos_from_start_and_strand(start, strand):
    try:
        s = int(float(start))
    except Exception:
        return None
    st = str(strand).strip()
    if st in ["+", "1", "1.0"]:
        return s + 16
    return s - 17


def _annotate_guides_with_output_region(guides_df, region_label="", region_start=None, region_end=None):
    if guides_df is None or guides_df.empty:
        return guides_df
    out = guides_df.copy()
    out["target_region_label"] = str(region_label) if region_label is not None else ""
    out["target_region_start"] = region_start if region_start is not None else ""
    out["target_region_end"] = region_end if region_end is not None else ""
    out["cut_pos"] = out.apply(lambda r: _cut_pos_from_start_and_strand(r.get("start", 0), r.get("strand", "+")), axis=1)
    if region_start is not None and region_end is not None:
        lo = min(int(region_start), int(region_end))
        hi = max(int(region_start), int(region_end))
        cut_series = pd.to_numeric(out["cut_pos"], errors="coerce")
        out["cut_in_target_region"] = cut_series.notna() & (cut_series >= lo) & (cut_series <= hi)
    else:
        out["cut_in_target_region"] = ""
    return out


def _safe_float(v, default=0.0):
    try:
        out = float(v)
        if pd.isna(out):
            return float(default)
        return out
    except Exception:
        return float(default)


def _get_loc2pos_types(chr_name, enst_id, pos, loc2posType):
    if loc2posType is None or len(loc2posType) == 0:
        return []
    chr_candidates = [str(chr_name)]
    if str(chr_name).lower().startswith("chr"):
        chr_candidates.append(str(chr_name)[3:])
    else:
        chr_candidates.append("chr" + str(chr_name))
    for cc in chr_candidates:
        if cc not in loc2posType:
            continue
        chr_dict = loc2posType[cc]
        if enst_id not in chr_dict:
            continue
        out = []
        mapping_dict = chr_dict[enst_id]
        for key in mapping_dict.keys():
            try:
                st, en = int(key[0]), int(key[1])
            except Exception:
                continue
            if int(pos) >= st and int(pos) <= en:
                out.append(mapping_dict[key])
        return out
    return []


def _psj_specificity_weight(specificity_score, low=45, high=65):
    x = _safe_float(specificity_score, 0.0)
    if x < low:
        return 0.0
    if x >= high:
        return 1.0
    return (x - low) / float(high - low)


def _psj_dist_weight(hdr_dist, variance=55):
    d = abs(_safe_float(hdr_dist, 0.0))
    return float(pow(2.718281828459045, (-1.0 * d * d) / (2.0 * variance)))


def _psj_position_weight(types):
    mapping = {
        "5UTR": 0.4,
        "3UTR": 1.0,
        "cds": 1.0,
        "within_2bp_of_exon_intron_junction": 0.01,
        "within_2bp_of_intron_exon_junction": 0.01,
        "3N4bp_up_of_exon_intron_junction": 0.1,
        "3_to_6bp_down_of_exon_intron_junction": 0.1,
        "3N4bp_up_of_intron_exon_junction": 0.1,
        "3N4bp_down_of_intron_exon_junction": 0.5,
    }
    if types is None or len(types) == 0:
        return 1.0
    lowest = 1.0
    for t in types:
        if t in mapping:
            lowest = min(lowest, float(mapping[t]))
    return lowest


def rerank_guides_csv_with_psj(
    guides_df,
    enst_id,
    chrom,
    insert_pos,
    enst_strand,
    terminus_type,
    loc2posType,
    spec_score_flavor,
    reg_penalty,
    alphas,
    max_abs_cut2ins=None,
    specificity_backend="default",
    chopchop_proxy_k=600.0,
    chopchop_proxy_n=1.0,
    preserve_input_order=False,
):
    if guides_df is None or guides_df.empty:
        return guides_df
    out = guides_df.copy().reset_index(drop=True)
    # Preserve raw CHOPCHOP ranking by default; keep PSJ reranking logic available behind non-default backends.
    if preserve_input_order and str(specificity_backend).lower() == "default":
        cut2ins = []
        for _, row in out.iterrows():
            cut_pos = _cut_pos_from_start_and_strand(row.get("start", 0), row.get("strand", "+"))
            if cut_pos is None:
                this_dist = 10**9
            else:
                this_dist = int(cut_pos) - int(insert_pos)
                if terminus_type == "start" and int(enst_strand) == -1:
                    this_dist += 1
                if terminus_type == "stop" and int(enst_strand) == 1:
                    this_dist += 1
            cut2ins.append(this_dist)
        out["chr"] = str(chrom)
        out["ID"] = str(enst_id)
        out["Insert_pos"] = int(insert_pos)
        out["Cut2Ins_dist"] = cut2ins
        out["spec_weight"] = float("nan")
        out["dist_weight"] = float("nan")
        out["pos_weight"] = float("nan")
        out["final_weight"] = float("nan")
        if max_abs_cut2ins is not None:
            out = out[pd.to_numeric(out["Cut2Ins_dist"], errors="coerce").abs() <= float(max_abs_cut2ins)]
        return out.reset_index(drop=True)
    spec_w = []
    dist_w = []
    pos_w = []
    final_w = []
    cut2ins = []
    for _, row in out.iterrows():
        cut_pos = _cut_pos_from_start_and_strand(row.get("start", 0), row.get("strand", "+"))
        if cut_pos is None:
            this_dist = 10**9
        else:
            this_dist = int(cut_pos) - int(insert_pos)
            if terminus_type == "start" and int(enst_strand) == -1:
                this_dist += 1
            if terminus_type == "stop" and int(enst_strand) == 1:
                this_dist += 1
        cut2ins.append(this_dist)

        if str(specificity_backend).lower() == "chopchop_proxy":
            penalty = _safe_float(row.get("chopchop_offtarget_penalty", float("nan")), float("nan"))
            if pd.isna(penalty):
                sw = 1.0
            else:
                sw = 1.0 / (1.0 + (float(penalty) / float(chopchop_proxy_k)) ** float(chopchop_proxy_n))
        else:
            sscore = _safe_float(row.get(spec_score_flavor, row.get("guideMITScore", 50.0)), 50.0)
            sw = _psj_specificity_weight(sscore)
        dw = _psj_dist_weight(this_dist)
        if reg_penalty:
            types_here = _get_loc2pos_types(chrom, enst_id, cut_pos, loc2posType) if cut_pos is not None else []
            types_next = _get_loc2pos_types(chrom, enst_id, int(cut_pos) + 1, loc2posType) if cut_pos is not None else []
            pw = min(_psj_position_weight(types_here), _psj_position_weight(types_next))
        else:
            pw = 1.0
        fscore = 1.0
        for w, a in zip([sw, dw, pw], alphas):
            if float(a) != 0:
                fscore *= float(w) ** float(a)
        spec_w.append(sw)
        dist_w.append(dw)
        pos_w.append(pw)
        final_w.append(fscore)

    out["chr"] = str(chrom)
    out["ID"] = str(enst_id)
    out["Insert_pos"] = int(insert_pos)
    out["Cut2Ins_dist"] = cut2ins
    out["spec_weight"] = spec_w
    out["dist_weight"] = dist_w
    out["pos_weight"] = pos_w
    out["final_weight"] = final_w

    if max_abs_cut2ins is not None:
        out = out[pd.to_numeric(out["Cut2Ins_dist"], errors="coerce").abs() <= float(max_abs_cut2ins)]

    if preserve_input_order:
        out = out.reset_index(drop=True)
    else:
        out = out.sort_values(["final_weight"], ascending=False, kind="stable").reset_index(drop=True)
    return out


def select_guides_from_csv(
    guides_df,
    entry,
    enst_id,
    terminus,
    chrom=None,
    coordinate=None,
    fallback_insert_pos=None,
    fallback_chr=None,
    max_abs_cut2ins=None,
):
    if guides_df is None or guides_df.empty:
        return pd.DataFrame()

    entry = str(entry)
    enst_id = str(enst_id)
    term = str(terminus)

    subset = pd.DataFrame()
    if "Entry" in guides_df.columns and "terminus" in guides_df.columns:
        subset = guides_df[
            (guides_df["Entry"].astype(str) == entry)
            & (guides_df["terminus"].astype(str) == term)
        ]
    if subset.empty and "ID" in guides_df.columns and "terminus" in guides_df.columns:
        subset = guides_df[
            (guides_df["ID"].astype(str) == enst_id)
            & (guides_df["terminus"].astype(str) == term)
        ]
    if subset.empty and term == "-" and chrom is not None and coordinate is not None and "chr" in guides_df.columns:
        subset = guides_df[guides_df["chr"].astype(str) == str(chrom)]

    subset = subset.copy().reset_index(drop=True)

    # Fallback: if ENST-mode N/C rows are absent, derive them from "-" rows.
    if (
        subset.empty
        and term in ["N", "C"]
        and (fallback_insert_pos is not None)
        and "terminus" in guides_df.columns
    ):
        if "ID" in guides_df.columns:
            base = guides_df[
                (guides_df["ID"].astype(str) == enst_id)
                & (guides_df["terminus"].astype(str).str.strip() == "-")
            ].copy()
        else:
            base = pd.DataFrame()
        if base.empty and "Entry" in guides_df.columns:
            base = guides_df[
                (guides_df["Entry"].astype(str) == entry)
                & (guides_df["terminus"].astype(str).str.strip() == "-")
            ].copy()
        if not base.empty:
            base["terminus"] = term
            base["Insert_pos"] = int(fallback_insert_pos)
            if fallback_chr is not None and "chr" in base.columns:
                base["chr"] = str(fallback_chr)
            cut2ins = []
            for _, r in base.iterrows():
                cp = _cut_pos_from_start_and_strand(r.get("start", 0), r.get("strand", "+"))
                if cp is None:
                    cut2ins.append(0)
                else:
                    cut2ins.append(int(cp) - int(fallback_insert_pos))
            base["Cut2Ins_dist"] = cut2ins
            subset = base.reset_index(drop=True)

    if (not subset.empty) and (max_abs_cut2ins is not None) and ("Cut2Ins_dist" in subset.columns):
        tmp = subset.copy()
        tmp["_abs_cut2ins"] = pd.to_numeric(tmp["Cut2Ins_dist"], errors="coerce").abs()
        tmp = tmp[tmp["_abs_cut2ins"].notna()]
        tmp = tmp[tmp["_abs_cut2ins"] <= float(max_abs_cut2ins)]
        subset = tmp.drop(columns=["_abs_cut2ins"], errors="ignore").reset_index(drop=True)

    return subset


def write_guide_table_rows(handle, entry, enst_id, terminus, ranked_df):
    if ranked_df is None or ranked_df.empty:
        return
    for rank, (_, row) in enumerate(ranked_df.iterrows(), start=1):
        handle.write(
            f"{entry},{enst_id},{terminus},{rank},{row.get('chopchop_rank','')},{row.get('chr','')},{row.get('Insert_pos','')},{row.get('seq','')},{row.get('pam','')},"
            f"{row.get('start','')},{row.get('end','')},{row.get('strand','')},{row.get('cut_pos','')},{row.get('target_region_label','')},{row.get('target_region_start','')},{row.get('target_region_end','')},"
            f"{row.get('Eff_scores','')},{row.get('MM0','')},{row.get('MM1','')},{row.get('MM2','')},{row.get('MM3','')},{row.get('Cut2Ins_dist','')}\n"
        )
    handle.flush()

def test_memory(n):
    """
    try allocate n MB of memory
    return true if can, and false otherwise
    """
    try:
        x = bytearray(1024 * 1000 * n)
        del x
        return True
    except:
        return False

def get_current_date_formatted():
    # Get the current date
    current_date = datetime.datetime.now()
    # Format the date as yyyy-mm-dd
    formatted_date = current_date.strftime("%d-%b-%Y").upper()
    return formatted_date

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print(
        'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(
            filename, lineno, line.strip(), exc_obj
        )
    )


if __name__ == "__main__":
    main()
