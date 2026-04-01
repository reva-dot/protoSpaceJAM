from Bio.Seq import Seq
import os
import os.path
import pandas as pd
from Bio.Seq import reverse_complement
import argparse
import sys
import math
import pickle
import logging
import time
import json
import subprocess
import tempfile
import re
import csv
from io import StringIO
from urllib.parse import urljoin
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

from protoSpaceJAM.util.hdr import HDR_flank #uncomment this for pip installation
# from util.hdr import HDR_flank

GUIDE_COLUMNS = [
    "chr",
    "seq",
    "pam",
    "start",
    "end",
    "strand",
    "guideMITScore",
    "guideCfdScore",
    "guideCfdScorev2",
    "guideCfdScorev3",
    "Eff_scores",
    "MM0",
    "MM1",
    "MM2",
    "MM3",
    "chopchop_offtarget_penalty",
    "chopchop_rank",
]

GUIDE_EXPORT_COLUMNS = GUIDE_COLUMNS + ["score_note"]

GENOME_TO_ENSEMBL = {
    "GRCh38": ("homo_sapiens", "GRCh38"),
    "GRCm39": ("mus_musculus", "GRCm39"),
    "GRCz11": ("danio_rerio", "GRCz11"),
    "mRatBN7.2": ("rattus_norvegicus", "mRatBN7.2"),
}

GENOME_TO_CHOPCHOP = {
    "GRCh38": "hg38",
    "GRCm39": "mm39",
    "GRCz11": "danRer11",
    "mRatBN7.2": "mRatBN7.2",
}


def _dump_chopchop_debug(
    chopchop_config,
    source,
    chrom,
    pos,
    window_start,
    window_end,
    target_seq,
    payload=None,
):
    if chopchop_config is None or not chopchop_config.get("debug_dump", False):
        return
    outdir = chopchop_config.get("debug_outdir", ".")
    debug_dir = os.path.join(outdir, "chopchop_debug")
    os.makedirs(debug_dir, exist_ok=True)
    stem = f"{source}_{str(chrom).replace(':', '_')}_{int(pos)}_{int(window_start)}_{int(window_end)}"

    fasta_path = os.path.join(debug_dir, f"{stem}.fa")
    with open(fasta_path, "w") as fh:
        fh.write(f">{chrom}:{window_start}-{window_end}\n{target_seq}\n")

    meta = {
        "source": source,
        "chrom": str(chrom),
        "pos": int(pos),
        "window_start": int(window_start),
        "window_end": int(window_end),
        "sequence_length": len(target_seq),
    }
    if payload is not None:
        meta["payload"] = payload
    with open(os.path.join(debug_dir, f"{stem}.json"), "w") as jh:
        json.dump(meta, jh, indent=2)

    # Also keep a flat CSV log for quick cross-checking against CHOPCHOP web runs.
    csv_path = os.path.join(debug_dir, "chopchop_debug_windows.csv")
    write_header = not os.path.isfile(csv_path)
    with open(csv_path, "a", newline="") as ch:
        writer = csv.writer(ch)
        if write_header:
            writer.writerow(
                [
                    "source",
                    "chrom",
                    "pos",
                    "window_start",
                    "window_end",
                    "sequence_length",
                    "sequence_window",
                    "fasta_path",
                    "json_path",
                ]
            )
        writer.writerow(
            [
                source,
                str(chrom),
                int(pos),
                int(window_start),
                int(window_end),
                len(target_seq),
                target_seq,
                fasta_path,
                os.path.join(debug_dir, f"{stem}.json"),
            ]
        )

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

#################
# custom logging #
#################

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
# The background is set with 40 plus the number of the color, and the foreground with 30

# These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"


def formatter_message(message, use_color=True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message


COLORS = {
    "WARNING": YELLOW,
    "INFO": WHITE,
    "DEBUG": BLUE,
    "CRITICAL": YELLOW,
    "ERROR": RED,
}


class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = (
                COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            )
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s "  # ($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)

    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler(stream=sys.stdout)
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.INFO)


def get_phase_in_codon(Chr, Pos, ENST_ID, ENST_PhaseInCodon):
    """
    get the phase in codon for the current position and flanking 2 bp
    returns a dictionary:
    {0:[-1/-2/-3/1/2/3/0] 0 stands for not coding sequence (current position is not in a codon)
    -1:[-1/-2/-3/1/2/3/0]
    -2:[-1/-2/-3/1/2/3/0]
    +1:[-1/-2/-3/1/2/3/0]
    +2:[-1/-2/-3/1/2/3/0]}
    """
    mydict = {-2: 0, -1: 0, 0: 0, 1: 0, 2: 0}
    if Chr in ENST_PhaseInCodon.keys():
        Chr_dict = ENST_PhaseInCodon[Chr]
        if Pos in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos].keys():
                mydict[0] = Chr_dict[Pos][ENST_ID]
        if Pos + 1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos + 1].keys():
                mydict[+1] = Chr_dict[Pos + 1][ENST_ID]
        if Pos + 2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos + 2].keys():
                mydict[+2] = Chr_dict[Pos + 2][ENST_ID]
        if Pos - 1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos - 1].keys():
                mydict[-1] = Chr_dict[Pos - 1][ENST_ID]
        if Pos - 2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos - 2].keys():
                mydict[-2] = Chr_dict[Pos - 2][ENST_ID]
    return mydict


def get_phase_in_codon0(Chr, Pos, ENST_ID, ENST_PhaseInCodon):
    """
    get the phase in codon for the current position and flanking 2 bp
    returns an int:
    -1/-2/-3/1/2/3/0 0 stands for not coding sequence (current position is not in a codon)
    """
    myInt = 0
    if Chr in ENST_PhaseInCodon.keys():
        Chr_dict = ENST_PhaseInCodon[Chr]
        if Pos in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos].keys():
                myInt = Chr_dict[Pos][ENST_ID]
    return myInt


def get_range(start, end):  # TODO phase should be ENST specific
    """
    return a list of number, start to end, step size = 1, -1 (if start > end)
    """
    if start <= end:
        return list(range(start, end + 1, 1))
    else:
        return list(range(start, end - 1, -1))


def get_HDR_template(
    df,
    ENST_info,
    type,
    ENST_PhaseInCodon,
    HDR_arm_len,
    genome_ver,
    payload_type, # define insertion or SNP
    tag, # for insertion
    SNP_payload, # for SNP
    loc2posType,
    ssODN_max_size,
    recoding_args,
    Donor_type,
    Strand_choice,
    syn_check_args,
    coordinate_without_ENST,
):
    for index, row in df.iterrows():
        ENST_ID = row["ID"]
        ENST_genename = ENST_info[ENST_ID].name
        ENST_strand = ENST_info[ENST_ID].features[0].strand
        Chr = row["chr"]
        InsPos = row[
            "Insert_pos"
        ]  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G", or the letter before the payload insertion in genomics-coord intersion mode,  in SNP mode, it is the first letter of the SNP payload replacement
        gStrand = convert_strand(row["strand"])  # gRNA strand
        guide_lo, guide_hi = _normalize_guide_bounds(row["start"], row["end"])
        gStart = guide_lo if gStrand == 1 else guide_hi
        CutPos = get_cut_pos(row["start"], row["end"], gStrand)
        Cut2Ins_dist = row["Cut2Ins_dist"]

        ##########################
        # important debug info
        # print(f"ENST: {ENST_ID} ENST_strand: {ENST_strand} chr {Chr} type {type}-tagging InsPos {InsPos} gStrand {gStrand} CutPos {CutPos} gStart {gStart} ")
        ##########################
        # CutPos_phase = get_phase_in_codon(Chr=Chr, Pos=InsPos, ENST_PhaseInCodon=ENST_PhaseInCodon)
        # print(f"CutPos phase: {CutPos_phase}")170

        # get target_seq
        # target_seq = get_target_seq(Chr= Chr, InsPos = InsPos, gRNAstrand = gStrand, CutPos = CutPos , type = type, ENST_ID = ENST_ID, ENST_info = ENST_info)

        leftArm, rightArm, left_start, left_end, right_start, right_end = get_HDR_arms(
            loc=[Chr, InsPos, ENST_strand],
            half_len=HDR_arm_len,
            type=type,
            genome_ver=genome_ver,
        )  # start>end is possible
        left_Arm_Phases = [
            get_phase_in_codon0(
                Chr=Chr, Pos=i, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon
            )
            for i in get_range(left_start, left_end)
        ]
        right_Arm_Phases = [
            get_phase_in_codon0(
                Chr=Chr, Pos=i, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon
            )
            for i in get_range(right_start, right_end)
        ]

        # support for coordinate without ENST
        if coordinate_without_ENST:
            left_Arm_Phases = [0 for i in left_Arm_Phases] # treat as non-coding
            right_Arm_Phases = [0 for i in right_Arm_Phases] # treat as non-coding

        myflank = HDR_flank(
            left_flk_seq=leftArm,
            right_flk_seq=rightArm,
            left_flk_coord_lst=[left_start, left_end],
            right_flk_coord_lst=[right_start, right_end],
            left_flk_phases=left_Arm_Phases,
            right_flk_phases=right_Arm_Phases,
            type=type,
            ENST_ID=ENST_ID,
            name=ENST_genename,
            ENST_strand=ENST_strand,
            ENST_chr=Chr,
            gStart=gStart,
            gStrand=gStrand,
            InsPos=InsPos,
            CutPos=CutPos,
            Cut2Ins_dist=Cut2Ins_dist,
            payload_type=payload_type,
            tag=tag,
            SNP_payload=SNP_payload,
            loc2posType=loc2posType,
            ssODN_max_size=ssODN_max_size,
            Donor_type=Donor_type,
            Strand_choice=Strand_choice,
            recoding_args=recoding_args,

            syn_check_args=syn_check_args,
            coordinate_without_ENST=coordinate_without_ENST,
        )
        return myflank
        # log IDs whose gRNA is not in the default-size HDR arms
        # if myflank.entire_gRNA_in_HDR_arms == False:
        #    gRNA_out_of_arms[type][ENST_ID] = False


# def get_target_seq(Chr, InsPos, gRNAstrand, CutPos, type, ENST_ID, ENST_info):
#     """
#     The target sequence is used to instantiate the HDR class (gdingle)
#     The target sequence should be in the direction of the gene. Reading from
#     left to right, it should have either a ATG or one of TAG, TGA, or TAA.
#     """
#     ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
#     if type == "start":
#         target_codon_loc = ATG_loc
#         target_codon_strand = ATG_loc[3]
#         if target_codon_strand == 1:
#             ATG_loc[2]
#     elif type == "stop":
#         target_codon_loc = stop_loc
#         target_codon_strand = stop_loc[3]
#     else:
#         sys.exit(f"unknown type: {type}")


def convert_strand(strand):
    # convert strand from +/- to 1/-1
    if strand == "+":
        return 1
    elif strand == "-":
        return -1
    else:
        return f"input strand:{strand} needs to be +/-"


def get_HDR_arms(loc, half_len, type, genome_ver):
    """
    input:  loc         [chr,pos,strand]  #start < end , strand is the coding strand
            half_len      length of the HDR arm (one sided)
    return: HDR arms -> in coding strand <-
            [5'arm, 3'arm]
    """
    Chr, Pos, Strand = loc
    # get arms
    vanilla_left_arm, vanilla_right_arm = None, None
    if type == "start":
        if Strand == 1:
            vanilla_left_arm = get_seq(
                chr=Chr,
                start=Pos - half_len + 1,
                end=Pos + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            vanilla_right_arm = get_seq(
                chr=Chr,
                start=Pos + 1,
                end=Pos + half_len + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            return [
                vanilla_left_arm,
                vanilla_right_arm,
                Pos - half_len + 1,
                Pos,
                Pos + 1,
                Pos + half_len,
            ]
        elif Strand == -1:
            vanilla_left_arm = get_seq(
                chr=Chr, start=Pos - half_len, end=Pos, strand=1, genome_ver=genome_ver
            )
            vanilla_right_arm = get_seq(
                chr=Chr, start=Pos, end=Pos + half_len, strand=1, genome_ver=genome_ver
            )
            return [
                reverse_complement(vanilla_right_arm),
                reverse_complement(vanilla_left_arm),
                Pos + half_len - 1,
                Pos,
                Pos - 1,
                Pos - half_len,
            ]
        else:
            sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")
    elif type == "stop":
        if Strand == 1:
            vanilla_left_arm = get_seq(
                chr=Chr, start=Pos - half_len, end=Pos, strand=1, genome_ver=genome_ver
            )
            vanilla_right_arm = get_seq(
                chr=Chr, start=Pos, end=Pos + half_len, strand=1, genome_ver=genome_ver
            )
            return [
                vanilla_left_arm,
                vanilla_right_arm,
                Pos - half_len,
                Pos - 1,
                Pos,
                Pos + half_len - 1,
            ]
        elif Strand == -1:
            vanilla_left_arm = get_seq(
                chr=Chr,
                start=Pos - half_len + 1,
                end=Pos + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            vanilla_right_arm = get_seq(
                chr=Chr,
                start=Pos + 1,
                end=Pos + half_len + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            return [
                reverse_complement(vanilla_right_arm),
                reverse_complement(vanilla_left_arm),
                Pos + half_len,
                Pos + 1,
                Pos,
                Pos - half_len + 1,
            ]
    else:
        sys.exit("unknown type {type}, acceptable values: start, stop")

    # rev_com ajustment
    if Strand == 1:
        return [vanilla_left_arm, vanilla_right_arm]
    elif Strand == -1:
        return [
            reverse_complement(vanilla_right_arm),
            reverse_complement(vanilla_left_arm),
        ]
    else:
        sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")

def get_gRNAs_target_coordinate(
    ENST_ID,
    chrom,
    pos,
    ENST_info,
    freq_dict,
    loc2file_index,
    loc2posType,
    genome_ver,
    pam,
    spec_score_flavor,
    reg_penalty,
    alphas,
    guide_source="precomputed",
    chopchop_config=None,
    dist=50,
):
    """
    input
        ENST_ID: ENST ID
        chr: chromosome (as part of the target location)
        pos: position (as part of the target location)
        ENST_info: gene model info loaded from pickle file
        loc2file_index: mapping of location to the file part that stores gRNAs in that region
        loc2posType: mapping of location to location type (e.g. 5UTR etc)
        dist: max cut to insert distance (default 50)
        reg_penalty: whether to penalize gRNAs that cut in the region of interest (default True)
    return:
        a dictionary of guide RNAs which cuts <[dist] to the target location

    """
    # get ENST strand
    ENST_strand = ENST_info[ENST_ID].features[0].strand
    ##################################
    # get gRNAs around the coordinate#
    ##################################
    # get the start codon chromosomal location
    target_pos = [chrom, pos, ENST_strand]

    # get gRNAs around the target location/coordinate
    df_gRNAs_target_pos = get_gRNAs_near_loc(
        loc=target_pos,
        dist=dist,
        loc2file_index=loc2file_index,
        genome_ver=genome_ver,
        pam=pam,
        guide_source=guide_source,
        chopchop_config=chopchop_config,
    )
    keep_chopchop_order = (
        str(guide_source).lower() == "chopchop"
        and str(((chopchop_config or {}).get("scorer_config", {}) or {}).get("backend", "default")).lower() == "default"
        and not (chopchop_config or {}).get("psj_rank_chopchop", False)
    )
    if keep_chopchop_order:
        ranked_df_gRNAs_ATG = _decorate_guides_for_hdr_without_reranking(
            loc=target_pos,
            gRNA_df=df_gRNAs_target_pos,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="start",
        )
    else:
        ranked_df_gRNAs_target_pos = rank_gRNAs_for_tagging(
            loc=target_pos,
            gRNA_df=df_gRNAs_target_pos,
            loc2posType=loc2posType,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="start",
            spec_score_flavor=spec_score_flavor,
            reg_penalty=reg_penalty,
            alphas=alphas,
            specificity_backend=str(((chopchop_config or {}).get("scorer_config", {}) or {}).get("backend", "default")),
            chopchop_proxy_k=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_k", 600.0)),
            chopchop_proxy_n=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_n", 1.0)),
        )
        ranked_df_gRNAs_ATG = ranked_df_gRNAs_target_pos.sort_values(
            "final_weight", ascending=False
        )  # sort descending on final weight

    return ranked_df_gRNAs_ATG


def get_gRNAs(
    ENST_ID,
    ENST_info,
    freq_dict,
    loc2file_index,
    loc2posType,
    genome_ver,
    pam,
    spec_score_flavor,
    reg_penalty,
    alphas,
    guide_source="precomputed",
    chopchop_config=None,
    dist=50,
):
    """
    Rank gRNAs that cut near the start and stop codon of an ENST_ID
    input
        ENST_ID: ENST ID
        ENST_info: gene model info loaded from pickle file
        loc2file_index: mapping of location to the file part that stores gRNAs in that region
        loc2posType: mapping of location to location type (e.g. 5UTR etc)
        dist: max cut to insert distance (default 50)
    return:
        a dictionary of guide RNAs which cuts <[dist] to end of start codon
        a dictionary of guide RNAs which cuts <[dist] to start of stop codon
    """
    # get location of start and stop location
    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
    # get ENST strand
    ENST_strand = ENST_info[ENST_ID].features[0].strand
    ##################################
    # get gRNAs around the start codon#
    ##################################
    # get the start codon chromosomal location
    log.debug(f"ATG_loc: {ATG_loc}")
    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc)  # [chr, pos,strand]
    log.debug(f"end of the ATG: {end_of_ATG_loc}")
    # get gRNA around the chromosomeal location (near ATG)
    df_gRNAs_ATG = get_gRNAs_near_loc(
        loc=end_of_ATG_loc,
        dist=dist,
        loc2file_index=loc2file_index,
        genome_ver=genome_ver,
        pam=pam,
        guide_source=guide_source,
        chopchop_config=chopchop_config,
    )
    start_interval, start_region_label = _terminal_cds_region_for_terminus(ENST_ID, ENST_info, "start")
    df_gRNAs_ATG = _filter_guides_to_terminal_cds(
        ENST_ID=ENST_ID,
        ENST_info=ENST_info,
        gRNA_df=df_gRNAs_ATG,
        terminus_type="start",
    )
    df_gRNAs_ATG = _annotate_guides_with_target_interval(df_gRNAs_ATG, start_region_label, start_interval)
    keep_chopchop_order = (
        str(guide_source).lower() == "chopchop"
        and str(((chopchop_config or {}).get("scorer_config", {}) or {}).get("backend", "default")).lower() == "default"
        and not (chopchop_config or {}).get("psj_rank_chopchop", False)
    )
    if keep_chopchop_order:
        ranked_df_gRNAs_ATG = _decorate_guides_for_hdr_without_reranking(
            loc=end_of_ATG_loc,
            gRNA_df=df_gRNAs_ATG,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="start",
        )
    else:
        ranked_df_gRNAs_ATG = rank_gRNAs_for_tagging(
            loc=end_of_ATG_loc,
            gRNA_df=df_gRNAs_ATG,
            loc2posType=loc2posType,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="start",
            spec_score_flavor=spec_score_flavor,
            reg_penalty=reg_penalty,
            alphas=alphas,
            specificity_backend=str(((chopchop_config or {}).get("scorer_config", {}) or {}).get("backend", "default")),
            chopchop_proxy_k=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_k", 600.0)),
            chopchop_proxy_n=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_n", 1.0)),
        )
        ranked_df_gRNAs_ATG = ranked_df_gRNAs_ATG.sort_values(
            "final_weight", ascending=False
        )  # sort descending on final weight

    ##################################
    # get gRNAs around the stop  codon#
    ##################################
    # get gRNAs around the stop codon
    log.debug(f"stop_loc: {stop_loc}")
    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
    log.debug(f"start of stop: {start_of_stop_loc}")  # [chr, pos,strand]
    # get gRNA around the chromosomeal location (near stop location)
    df_gRNAs_stop = get_gRNAs_near_loc(
        loc=start_of_stop_loc,
        dist=dist,
        loc2file_index=loc2file_index,
        genome_ver=genome_ver,
        pam=pam,
        guide_source=guide_source,
        chopchop_config=chopchop_config,
    )
    stop_interval, stop_region_label = _terminal_cds_region_for_terminus(ENST_ID, ENST_info, "stop")
    df_gRNAs_stop = _filter_guides_to_terminal_cds(
        ENST_ID=ENST_ID,
        ENST_info=ENST_info,
        gRNA_df=df_gRNAs_stop,
        terminus_type="stop",
    )
    df_gRNAs_stop = _annotate_guides_with_target_interval(df_gRNAs_stop, stop_region_label, stop_interval)
    if keep_chopchop_order:
        ranked_df_gRNAs_stop = _decorate_guides_for_hdr_without_reranking(
            loc=start_of_stop_loc,
            gRNA_df=df_gRNAs_stop,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="stop",
        )
    else:
        ranked_df_gRNAs_stop = rank_gRNAs_for_tagging(
            loc=start_of_stop_loc,
            gRNA_df=df_gRNAs_stop,
            loc2posType=loc2posType,
            ENST_ID=ENST_ID,
            ENST_strand=ENST_strand,
            type="stop",
            spec_score_flavor=spec_score_flavor,
            reg_penalty=reg_penalty,
            alphas=alphas,
            specificity_backend=str(((chopchop_config or {}).get("scorer_config", {}) or {}).get("backend", "default")),
            chopchop_proxy_k=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_k", 600.0)),
            chopchop_proxy_n=float(((chopchop_config or {}).get("scorer_config", {}) or {}).get("chopchop_proxy_n", 1.0)),
        )
        ranked_df_gRNAs_stop = ranked_df_gRNAs_stop.sort_values(
            "final_weight", ascending=False
        )  # sort descending on final weight

    return [ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop]


def rank_gRNAs_for_tagging(
    loc, gRNA_df, loc2posType, ENST_ID, ENST_strand, type, spec_score_flavor, reg_penalty, alphas=[1, 1, 1], specificity_backend="default", chopchop_proxy_k=600.0, chopchop_proxy_n=1.0
):
    """
    input:  loc         [chr,pos,strand]  #start < end
            gRNA_df     pandas dataframe, *unranked*   columns: "seq","pam","start","end", "strand", "guideMITScore","guideCfdScore","guideCfdScorev2","guideCfdScorev3", "Eff_scores"  !! neg strand: start > end
            alphas       scaling factor for the weights, default [1,1,1]
            type        "start" or "stop:
    output: gRNA_df     pandas dataframe *ranked*      columns: "seq","pam","start","end", "strand", "guideMITScore","guideCfdScore","guideCfdScorev2","guideCfdScorev3", "Eff_scores"  !! neg strand: start > end
    """
    insPos = loc[
        1
    ]  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
    Chr = loc[0]

    col_spec_weight = []
    col_dist_weight = []
    col_pos_weight = []
    col_final_weight = []
    Chrs = []
    ENSTs = []
    InsertPos = []
    cut2insDist_list = []
    # assign a score to each gRNA
    for index, row in gRNA_df.iterrows():
        start = row[2]
        end = row[3]
        strand = row[4]
        # Get cut to insert distance
        cutPos = get_cut_pos(start, end, strand)
        cut2insDist = cutPos - insPos

        # adjust cut2insDist
        if type == "start" and ENST_strand == -1:
            cut2insDist += 1
        if type == "stop" and ENST_strand == 1:
            cut2insDist += 1

        # calc. specificity_weight
        CSS = row[spec_score_flavor]
        if str(specificity_backend).lower() == "chopchop_proxy":
            penalty = pd.to_numeric(row.get("chopchop_offtarget_penalty", float("nan")), errors="coerce")
            if pd.isna(penalty):
                specificity_weight = 1.0
            else:
                specificity_weight = 1.0 / (1.0 + (float(penalty) / float(chopchop_proxy_k)) ** float(chopchop_proxy_n))
        else:
            specificity_weight = _specificity_weight(CSS)
        col_spec_weight.append(specificity_weight)

        # calc. distance_weight
        distance_weight = _dist_weight(hdr_dist=cut2insDist)
        col_dist_weight.append(distance_weight)

        # get position_weight
        position_type = _get_position_type(
                chr=Chr, ID=ENST_ID, pos=cutPos, loc2posType=loc2posType
            )
        if reg_penalty is True:
            position_weight = _position_weight(position_type)

            position_type_nextbp = _get_position_type(
                chr=Chr, ID=ENST_ID, pos=cutPos + 1, loc2posType=loc2posType
            )
            position_weight_nextbp = _position_weight(position_type_nextbp)

            position_weight = min([position_weight, position_weight_nextbp])
        else:
            position_weight = 1

        col_pos_weight.append(position_weight)

        # add info to the df
        Chrs.append(Chr)
        ENSTs.append(ENST_ID)
        InsertPos.append(insPos)
        cut2insDist_list.append(cut2insDist)

        log.debug(
            f"strand {strand} {start}-{end} cutPos {cutPos} insert_loc {loc} cut2insDist {cut2insDist} distance_weight {distance_weight:.2f} CFD_score {CSS} specificity_weight {specificity_weight} pos_type {position_type} position_weight {position_weight}"
        )

        # calc. final_weight
        final_score = 1
        for weight, alpha in zip([specificity_weight, distance_weight, position_weight], alphas):
            if alpha != 0:
                final_score *= float(weight ** alpha)

        col_final_weight.append(final_score)

    # add info and weight columns to the df
    gRNA_df["chr"] = Chrs
    gRNA_df["ID"] = ENSTs
    gRNA_df["Insert_pos"] = InsertPos
    gRNA_df["Cut2Ins_dist"] = cut2insDist_list
    gRNA_df["spec_weight"] = col_spec_weight
    gRNA_df["dist_weight"] = col_dist_weight
    gRNA_df["pos_weight"] = col_pos_weight
    gRNA_df["final_weight"] = col_final_weight

    # rank gRNAs based on the score
    gRNA_df["final_pct_rank"] = gRNA_df["final_weight"].rank(pct=True)

    return gRNA_df


def _decorate_guides_for_hdr_without_reranking(loc, gRNA_df, ENST_ID, ENST_strand, type):
    """
    Keep incoming guide order (e.g. CHOPCHOP rank order), but add the
    columns required by downstream HDR generation and reporting.
    """
    if gRNA_df is None or gRNA_df.empty:
        return gRNA_df

    insPos = loc[1]
    Chr = loc[0]
    gRNA_df = gRNA_df.copy()

    cut2insDist_list = []
    for _, row in gRNA_df.iterrows():
        cutPos = get_cut_pos(row["start"], row["end"], row["strand"])
        cut2insDist = cutPos - insPos
        if type == "start" and ENST_strand == -1:
            cut2insDist += 1
        if type == "stop" and ENST_strand == 1:
            cut2insDist += 1
        cut2insDist_list.append(cut2insDist)

    gRNA_df["chr"] = Chr
    gRNA_df["ID"] = ENST_ID
    gRNA_df["Insert_pos"] = insPos
    gRNA_df["Cut2Ins_dist"] = cut2insDist_list
    # Preserve CHOPCHOP order while marking PSJ-only ranking fields as not computed.
    gRNA_df["spec_weight"] = float("nan")
    gRNA_df["dist_weight"] = float("nan")
    gRNA_df["pos_weight"] = float("nan")
    gRNA_df["final_weight"] = float("nan")
    gRNA_df["final_pct_rank"] = float("nan")
    return gRNA_df


def _normalize_guide_bounds(start, end):
    try:
        s = int(float(start))
        e = int(float(end))
    except Exception:
        return (None, None)
    return (min(s, e), max(s, e))


def get_cut_pos(start, end, strand):
    """
    start/end:gRNA genomic span, normalized or legacy
    strand:gRNA strand
    Returns the genomic base immediately upstream of the cut using the
    protoSpaceJAM historical convention.
    """
    lo, hi = _normalize_guide_bounds(start, end)
    if lo is None or hi is None:
        return None
    if strand == "+" or strand == "1" or strand == 1:
        return lo + 16
    return hi - 17


def _terminal_cds_region_for_terminus(ENST_ID, ENST_info, terminus_type):
    my_transcript = ENST_info[ENST_ID]
    cds_list = [feat for feat in my_transcript.features if feat.type == "CDS"]
    if len(cds_list) == 0:
        return (None, "CDS")
    target_index = 0 if terminus_type == "start" else (len(cds_list) - 1)
    target_cds = cds_list[target_index]
    start = int(target_cds.location.start)
    end = int(target_cds.location.end)
    interval = [min(start, end), max(start, end)]
    label = "CDS_exon1" if terminus_type == "start" else "CDS_terminal_exon"
    return interval, label


def _filter_guides_to_terminal_cds(ENST_ID, ENST_info, gRNA_df, terminus_type):
    if gRNA_df is None or gRNA_df.empty:
        return gRNA_df
    interval, _ = _terminal_cds_region_for_terminus(ENST_ID, ENST_info, terminus_type)
    if interval is None:
        return gRNA_df
    lo, hi = int(interval[0]), int(interval[1])
    tmp = gRNA_df.copy()
    tmp["_cut_pos"] = tmp.apply(lambda r: get_cut_pos(r["start"], r["end"], r["strand"]), axis=1)
    tmp = tmp[
        pd.to_numeric(tmp["_cut_pos"], errors="coerce").notna()
        & (pd.to_numeric(tmp["_cut_pos"], errors="coerce") >= lo)
        & (pd.to_numeric(tmp["_cut_pos"], errors="coerce") <= hi)
    ]
    return tmp.drop(columns=["_cut_pos"], errors="ignore")


def _annotate_guides_with_target_interval(gRNA_df, region_label, interval):
    if gRNA_df is None or gRNA_df.empty:
        return gRNA_df
    out = gRNA_df.copy()
    out["target_region_label"] = str(region_label)
    if interval is None:
        out["target_region_start"] = float("nan")
        out["target_region_end"] = float("nan")
        out["cut_pos"] = out.apply(lambda r: get_cut_pos(r["start"], r["end"], r["strand"]), axis=1)
        out["cut_in_target_region"] = float("nan")
        return out
    lo, hi = int(interval[0]), int(interval[1])
    out["target_region_start"] = lo
    out["target_region_end"] = hi
    out["cut_pos"] = out.apply(lambda r: get_cut_pos(r["start"], r["end"], r["strand"]), axis=1)
    out["cut_in_target_region"] = (
        pd.to_numeric(out["cut_pos"], errors="coerce").notna()
        & (pd.to_numeric(out["cut_pos"], errors="coerce") >= lo)
        & (pd.to_numeric(out["cut_pos"], errors="coerce") <= hi)
    )
    return out


def _get_position_type(chr, ID, pos, loc2posType):
    """
    #input: mostly self-explanatory, loc2posType is a dictionary that translates location into types (e.g. exon/intron junctions etc)
    #return a list of types for the input position/ID combination
    >>> loc2posType = read_pickle_files(os.path.join("..","genome_files","parsed_gff3", "GRCh38","loc2posType.pickle"))

    #intron - exon junction
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399367, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399368, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['within_3bp_of_intron_exon_junction', '3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399369, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399370, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399371, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399372, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399373, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_intron_exon_junction', '3N4bp_down_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399374, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['cds', '3N4bp_down_of_intron_exon_junction']

    #intron - exon junction
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399481, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['cds', '3N4bp_up_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399482, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_exon_intron_junction', '3N4bp_up_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399483, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399484, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['cds', 'within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399485, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399486, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399487, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['within_3bp_of_exon_intron_junction', '3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399488, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399489, loc2posType = loc2posType) #5bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399490, loc2posType = loc2posType) #6bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']

    #non cds exon
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398847, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398848, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['within_3bp_of_intron_exon_junction', '3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398849, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398850, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398851, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['5UTR', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398852, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398853, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_intron_exon_junction', '3N4bp_down_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398854, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['cds', '3N4bp_down_of_intron_exon_junction']
    """
    # print(f"{type(chr)} {ID} {type(pos)}")
    if chr not in loc2posType:
        return []
    chr_dict = loc2posType[chr]
    if not ID in chr_dict.keys():
        return []
    else:
        types = []
        mapping_dict = chr_dict[ID]
        for key in mapping_dict.keys():
            if in_interval_leftrightInclusive(pos, key):
                types.append(mapping_dict[key])
        return types


def _position_weight(types):
    """
    input: a list of types
    output: the lowest weight among all the types
    """
    mapping_dict = {
        "5UTR": 0.4,
        "3UTR": 1,
        "cds": 1,
        "within_2bp_of_exon_intron_junction": 0.01,
        "within_2bp_of_intron_exon_junction": 0.01,
        "3N4bp_up_of_exon_intron_junction": 0.1,
        "3_to_6bp_down_of_exon_intron_junction": 0.1,
        "3N4bp_up_of_intron_exon_junction": 0.1,
        "3N4bp_down_of_intron_exon_junction": 0.5,
    }
    lowest_weight = 1
    for type in types:
        if type in mapping_dict.keys():
            weight = mapping_dict[type]
            if weight < lowest_weight:
                lowest_weight = weight
        else:
            pass
            # sys.exit(f"unexpected position type: {type}")
    return lowest_weight


def _dist_weight(hdr_dist: int, _dist_weight_variance=55) -> float:
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _dist_weight(0)
    1.0
    >>> _dist_weight(5)
    0.7967034698934616
    >>> _dist_weight(10)
    0.402890321529133
    >>> _dist_weight(-20)
    0.026347980814448734
    """
    variance = _dist_weight_variance

    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist ** 2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight


def _specificity_weight(
    specificity_score: float, _specificity_weight_low=45, _specificity_weight_high=65
):
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _specificity_weight(20)
    0
    >>> _specificity_weight(60)
    0.75
    >>> _specificity_weight(80)
    1
    """
    low = _specificity_weight_low
    high = _specificity_weight_high

    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= low:
        return 0
    elif specificity_score >= high:
        return 1
    else:
        return 1 / (high - low) * (specificity_score - low)


def get_end_pos_of_ATG(ATG_loc):
    """
    input: ATG_loc              [chr,start,end,strand] #start < end
    output: the pos of G in ATG [chr,pos,strand]       #start < end
    """
    strand = ATG_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [ATG_loc[0], ATG_loc[2], ATG_loc[3]]
    else:
        return [ATG_loc[0], ATG_loc[1], ATG_loc[3]]


def get_start_pos_of_stop(stop_loc):
    """
    input: stop_loc                                 [chr,start,end,strand]  #start < end
    output: the pos of first base in the stop codon [chr,pos,strand]        #start < end
    """
    strand = stop_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [stop_loc[0], stop_loc[1], stop_loc[3]]
    else:
        return [stop_loc[0], stop_loc[2], stop_loc[3]]


def get_start_stop_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
    output: a list of three items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
            Exon_end_ATG: Bool
    """
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
    # check if ATG is at the end of the first exon
    if check_ATG_at_exonEnd(my_transcript):
        CDS_first = cdsList[
            1
        ]  # use the second cds if ATG is at the end of the first exon
    # cdsList is in transcript order, so the first CDS contains the start codon
    # and the last CDS contains the stop codon for both strands.
    if CDS_first.strand == 1:
        ATG_loc = [
            CDS_first.location.ref,
            CDS_first.location.start + 0,
            CDS_first.location.start + 2,
            1,
        ]
    else:
        ATG_loc = [
            CDS_first.location.ref,
            CDS_first.location.end - 2,
            CDS_first.location.end + 0,
            -1,
        ]

    if CDS_last.strand == 1:
        stop_loc = [
            CDS_last.location.ref,
            CDS_last.location.end - 2,
            CDS_last.location.end + 0,
            1,
        ]
    else:
        stop_loc = [
            CDS_last.location.ref,
            CDS_last.location.start + 0,
            CDS_last.location.start + 2,
            -1,
        ]
    return [ATG_loc, stop_loc]


def _empty_guides_df():
    return pd.DataFrame(columns=GUIDE_COLUMNS)


def _canonize_chromosome(chr_name):
    chr_name = str(chr_name)
    if chr_name.lower().startswith("chr"):
        return chr_name[3:]
    return chr_name


def _build_ensembl_sequence_url(chrom, start, end, strand, genome_ver):
    species, _assembly = GENOME_TO_ENSEMBL.get(genome_ver, ("homo_sapiens", "GRCh38"))
    chrom = _canonize_chromosome(chrom)
    strand_num = 1 if str(strand) in ("1", "+", "plus") else -1
    region = f"{chrom}:{int(start)}..{int(end)}:{strand_num}"
    return f"https://rest.ensembl.org/sequence/region/{species}/{region}"


def fetch_sequence_from_ensembl(chrom, start, end, strand, genome_ver, timeout=30, max_retries=5):
    url = _build_ensembl_sequence_url(chrom, start, end, strand, genome_ver)
    last_exc = None
    for attempt in range(max(1, int(max_retries))):
        req = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "protoSpaceJAM/ensembl-seq-fetch",
            },
        )
        try:
            with urlopen(req, timeout=timeout) as response:
                data = json.loads(response.read().decode("utf-8"))
            if "seq" not in data:
                raise RuntimeError(
                    f"Ensembl sequence response missing 'seq' for {chrom}:{start}-{end}"
                )
            return data["seq"].upper()
        except (HTTPError, URLError, TimeoutError, ValueError, RuntimeError) as exc:
            last_exc = exc
            is_http_5xx = isinstance(exc, HTTPError) and 500 <= int(exc.code) < 600
            is_retryable = is_http_5xx or isinstance(exc, (URLError, TimeoutError))
            if attempt >= int(max_retries) - 1 or not is_retryable:
                break
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"Unable to fetch sequence from Ensembl ({url}): {last_exc}") from last_exc


def _best_numeric_series(df, candidates, default_value):
    for name in candidates:
        col_name = _find_column_case_insensitive(df, [name])
        if col_name is not None:
            out = pd.to_numeric(df[col_name], errors="coerce")
            if out.notna().any():
                return out
    return pd.Series([default_value] * len(df), index=df.index, dtype=float)


def _find_column_case_insensitive(df, candidate_names):
    cols = {str(c).lower().strip(): c for c in df.columns}
    for name in candidate_names:
        key = str(name).lower().strip()
        if key in cols:
            return cols[key]
    return None


def _extract_chopchop_offtarget_penalty(df_raw):
    mm0_col = _find_column_case_insensitive(df_raw, ["MM0"])
    mm1_col = _find_column_case_insensitive(df_raw, ["MM1"])
    mm2_col = _find_column_case_insensitive(df_raw, ["MM2"])
    mm3_col = _find_column_case_insensitive(df_raw, ["MM3"])
    if all(col is None for col in [mm0_col, mm1_col, mm2_col, mm3_col]):
        return pd.Series([float("nan")] * len(df_raw), index=df_raw.index, dtype=float)

    mm0 = pd.to_numeric(df_raw[mm0_col], errors="coerce") if mm0_col is not None else 0.0
    mm1 = pd.to_numeric(df_raw[mm1_col], errors="coerce") if mm1_col is not None else 0.0
    mm2 = pd.to_numeric(df_raw[mm2_col], errors="coerce") if mm2_col is not None else 0.0
    mm3 = pd.to_numeric(df_raw[mm3_col], errors="coerce") if mm3_col is not None else 0.0

    if not isinstance(mm0, pd.Series):
        mm0 = pd.Series([mm0] * len(df_raw), index=df_raw.index, dtype=float)
    if not isinstance(mm1, pd.Series):
        mm1 = pd.Series([mm1] * len(df_raw), index=df_raw.index, dtype=float)
    if not isinstance(mm2, pd.Series):
        mm2 = pd.Series([mm2] * len(df_raw), index=df_raw.index, dtype=float)
    if not isinstance(mm3, pd.Series):
        mm3 = pd.Series([mm3] * len(df_raw), index=df_raw.index, dtype=float)

    penalty = (
        1000.0 * mm0.fillna(0.0)
        + 800.0 * mm1.fillna(0.0)
        + 600.0 * mm2.fillna(0.0)
        + 400.0 * mm3.fillna(0.0)
    )
    return penalty.astype(float)


def _run_crispor_scores(df_guides, genome_ver, pam, scorer_config):
    if scorer_config is None or str(scorer_config.get("backend", "default")).lower() != "crispor":
        return None

    cmd_template = str(scorer_config.get("cmd_template", "")).strip()
    if cmd_template == "":
        raise RuntimeError(
            "specificity_backend=crispor requires --crispor_cmd_template."
        )

    df_input = df_guides.copy()
    df_input["seq"] = df_input["seq"].astype(str).str.upper()
    df_input["pam"] = df_input["pam"].astype(str).str.upper()
    max_guides = scorer_config.get("max_guides", None) if scorer_config is not None else None
    try:
        if max_guides is not None:
            max_guides = int(max_guides)
    except Exception:
        max_guides = None
    if max_guides is not None and max_guides > 0:
        df_input = df_input.head(max_guides).copy()

    with tempfile.TemporaryDirectory(prefix="protospacejam_crispor_") as temp_dir:
        input_file = os.path.join(temp_dir, "guides_for_crispor.tsv")
        output_file = os.path.join(temp_dir, "crispor_scores.tsv")
        export_cols = [c for c in ["chr", "seq", "pam", "start", "end", "strand"] if c in df_input.columns]
        df_input[export_cols].to_csv(
            input_file,
            sep="\t",
            index=False,
        )
        command = cmd_template.format(
            input_file=input_file,
            output_file=output_file,
            output_dir=temp_dir,
            genome=genome_ver or "",
            pam=pam,
        )
        proc = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            err_text = (proc.stderr or "").strip()
            out_text = (proc.stdout or "").strip()
            detail = "\n".join([x for x in [err_text, out_text] if x])
            raise RuntimeError(
                "CRISPOR scoring command failed with exit code "
                f"{proc.returncode}: {detail}"
            )
        if not os.path.isfile(output_file):
            raise RuntimeError(
                f"CRISPOR scoring command did not produce expected output file: {output_file}"
            )
        df_scores = pd.read_csv(output_file, sep=None, engine="python", dtype=str)

    if df_scores.empty:
        return None

    seq_col = _find_column_case_insensitive(
        df_scores,
        ["seq", "guide", "guideSeq", "target", "sgrna", "grna_seq"],
    )
    pam_col = _find_column_case_insensitive(
        df_scores,
        ["pam", "pamSeq", "pam_sequence"],
    )
    mit_col = _find_column_case_insensitive(
        df_scores,
        ["guideMITScore", "mitSpecScore", "mit_score", "mitscore", "specificity"],
    )
    cfd_col = _find_column_case_insensitive(
        df_scores,
        ["guideCfdScore", "cfdSpecScore", "cfd_score", "cfdscore"],
    )
    if seq_col is None:
        raise RuntimeError("CRISPOR score output must include a guide sequence column.")
    if mit_col is None and cfd_col is None:
        raise RuntimeError(
            "CRISPOR score output must include at least one MIT or CFD score column."
        )

    df_scores = df_scores.copy()
    df_scores["_seq_key"] = (
        df_scores[seq_col].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)
    )
    if pam_col is not None:
        df_scores["_pam_key"] = (
            df_scores[pam_col].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)
        )
    else:
        df_scores["_pam_key"] = ""

    keep_cols = ["_seq_key", "_pam_key"]
    if mit_col is not None:
        df_scores["guideMITScore"] = pd.to_numeric(df_scores[mit_col], errors="coerce")
        keep_cols.append("guideMITScore")
    if cfd_col is not None:
        df_scores["guideCfdScore"] = pd.to_numeric(df_scores[cfd_col], errors="coerce")
        keep_cols.append("guideCfdScore")
    return df_scores[keep_cols].drop_duplicates(subset=["_seq_key", "_pam_key"], keep="first")


def _apply_external_specificity_scores(df_guides, genome_ver, pam, scorer_config):
    df = df_guides.copy()
    mit_scores = _best_numeric_series(
        df,
        ["guideMITScore", "mitscore", "mit_score", "specificity", "offtargetscore"],
        float("nan"),
    )
    cfd_scores = _best_numeric_series(
        df,
        ["guideCfdScore", "cfdscore", "cfd_score", "cfdspecscore"],
        float("nan"),
    )
    score_note = ""

    backend = str((scorer_config or {}).get("backend", "default")).lower()

    if backend == "crispor":
        df["_seq_key"] = df["seq"].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)
        df["_pam_key"] = df["pam"].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)
        df_scores = _run_crispor_scores(df, genome_ver=genome_ver, pam=pam, scorer_config=scorer_config)
        if df_scores is not None and not df_scores.empty:
            merged = df.merge(df_scores, how="left", on=["_seq_key", "_pam_key"], suffixes=("", "_crispor"))
            if "guideMITScore_crispor" in merged.columns:
                mit_scores = pd.to_numeric(merged["guideMITScore_crispor"], errors="coerce")
            if "guideCfdScore_crispor" in merged.columns:
                cfd_scores = pd.to_numeric(merged["guideCfdScore_crispor"], errors="coerce")
            df = merged
            df = df.drop(columns=[c for c in ["guideMITScore_crispor", "guideCfdScore_crispor"] if c in df.columns])
            note_parts = []
            if mit_scores.notna().any():
                note_parts.append("guideMITScore=CRISPOR")
            if cfd_scores.notna().any():
                note_parts.append("guideCfdScore=CRISPOR")
            score_note = "; ".join(note_parts)
        df = df.drop(columns=[c for c in ["_seq_key", "_pam_key"] if c in df.columns])

    if backend == "chopchop_proxy":
        score_note = "spec_weight=CHOPCHOP_proxy_penalty; guideMITScore=compat_default"
    elif score_note == "":
        if mit_scores.notna().any():
            score_note = "guideMITScore=source_table"
        else:
            score_note = "guideMITScore=not_computed"

    df["guideMITScore"] = pd.to_numeric(mit_scores, errors="coerce")
    df["guideCfdScore"] = cfd_scores
    df["guideCfdScorev2"] = pd.Series([float("nan")] * len(df), index=df.index, dtype=float)
    df["guideCfdScorev3"] = pd.Series([float("nan")] * len(df), index=df.index, dtype=float)
    df["score_note"] = score_note
    return df


def _standardize_chopchop_df(df_raw, chr_name, pam, window_start):
    if df_raw is None or df_raw.empty:
        return _empty_guides_df()

    cols = {c.lower().strip(): c for c in df_raw.columns}
    seq_col = None
    for key in ["seq", "sequence", "guide", "target", "targetsequence", "target sequence", "sgrna"]:
        if key.lower() in cols:
            seq_col = cols[key.lower()]
            break
    if seq_col is None:
        raise RuntimeError("CHOPCHOP output is missing a guide sequence column")

    start_col = None
    end_col = None
    strand_col = None
    pam_col = None
    for key in ["start", "chromstart", "genomicstart", "from"]:
        if key.lower() in cols:
            start_col = cols[key.lower()]
            break
    for key in ["end", "chromend", "genomicend", "to"]:
        if key.lower() in cols:
            end_col = cols[key.lower()]
            break
    for key in ["strand", "orientation"]:
        if key.lower() in cols:
            strand_col = cols[key.lower()]
            break
    for key in ["pam", "pamseq", "pam_sequence"]:
        if key.lower() in cols:
            pam_col = cols[key.lower()]
            break

    rank_col = None
    for key in ["rank", "ranking"]:
        if key.lower() in cols:
            rank_col = cols[key.lower()]
            break

    genomic_loc_col = None
    for key in ["genomic location", "genomic_location", "location", "position"]:
        if key.lower() in cols:
            genomic_loc_col = cols[key.lower()]
            break

    if (start_col is None or end_col is None) and genomic_loc_col is None:
        raise RuntimeError("CHOPCHOP output must include start/end or genomic location columns")
    if strand_col is None:
        raise RuntimeError("CHOPCHOP output must include strand column")

    if start_col is not None and end_col is not None:
        start = pd.to_numeric(df_raw[start_col], errors="coerce")
        end = pd.to_numeric(df_raw[end_col], errors="coerce")
    else:
        # results.tsv exposes one genomic position; approximate start/end from that anchor.
        anchor = (
            df_raw[genomic_loc_col]
            .astype(str)
            .str.extract(r"(?P<chr>[^:]+):(?P<pos>\d+)")["pos"]
        )
        anchor = pd.to_numeric(anchor, errors="coerce")
        start = anchor.copy()
        end = anchor.copy()

    strand = df_raw[strand_col].astype(str).str.strip().replace({"1": "+", "-1": "-"})

    full_seq = df_raw[seq_col].astype(str).str.upper().str.replace(r"[^ACGTN]", "", regex=True)

    # CHOPCHOP outputs may be local window coordinates; convert when needed.
    max_coord = pd.concat([start, end], axis=1).max(axis=1)
    if max_coord.max(skipna=True) <= 5000:
        start = start + int(window_start) - 1
        end = end + int(window_start) - 1

    if start_col is None or end_col is None:
        # Approximate a genomic-low -> genomic-high protospacer span when only one
        # genomic location is provided.
        seq_len = full_seq.str.len().where(full_seq.str.len() > 0, 20)
        end = start + seq_len - 1

    df = pd.DataFrame(index=df_raw.index)
    df["chr"] = str(chr_name) if chr_name is not None else ""
    if genomic_loc_col is not None:
        extracted_chr = (
            df_raw[genomic_loc_col]
            .astype(str)
            .str.extract(r"(?P<chr>[^:]+):(?P<pos>\d+)")["chr"]
        )
        df["chr"] = extracted_chr.fillna(df["chr"])
    seq_guess = full_seq.where(full_seq.str.len() <= 20, full_seq.str.slice(0, 20))
    pam_guess = full_seq.where(full_seq.str.len() < 23, full_seq.str.slice(-3))
    df["seq"] = seq_guess
    df["pam"] = df_raw[pam_col].astype(str).str.upper() if pam_col else pam_guess.fillna(str(pam).upper())
    df["start"] = pd.concat([start, end], axis=1).min(axis=1).round().astype("Int64")
    df["end"] = pd.concat([start, end], axis=1).max(axis=1).round().astype("Int64")
    df["strand"] = strand

    eff_scores = _best_numeric_series(
        df_raw,
        ["eff_score", "efficiency", "doench", "eff_scores"],
        0.0,
    )
    df["Eff_scores"] = eff_scores
    for mm_name in ["MM0", "MM1", "MM2", "MM3"]:
        mm_col = _find_column_case_insensitive(df_raw, [mm_name])
        if mm_col is not None:
            df[mm_name] = pd.to_numeric(df_raw[mm_col], errors="coerce")
        else:
            df[mm_name] = pd.Series([float("nan")] * len(df), index=df.index, dtype=float)
    df["chopchop_offtarget_penalty"] = _extract_chopchop_offtarget_penalty(df_raw)
    if rank_col is not None:
        df["chopchop_rank"] = pd.to_numeric(df_raw[rank_col], errors="coerce")
    else:
        df["chopchop_rank"] = pd.Series(range(1, len(df) + 1), index=df.index, dtype=float)
    df = _apply_external_specificity_scores(
        df,
        genome_ver=None,
        pam=pam,
        scorer_config=None,
    )

    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["strand"] = df["strand"].where(df["strand"].isin(["+", "-"]), "+")
    return df[GUIDE_EXPORT_COLUMNS]


def convert_chopchop_raw_to_psj(
    df_raw,
    pam,
    window_start=1,
    desired_insert_pos=None,
    default_chr=None,
    genome_ver=None,
    scorer_config=None,
):
    """
    Convert raw CHOPCHOP results.tsv table to protoSpaceJAM guide schema.
    Keeps CHOPCHOP order.
    """
    df = _standardize_chopchop_df(
        df_raw=df_raw,
        chr_name=default_chr or "",
        pam=pam,
        window_start=window_start,
    )
    df = _apply_external_specificity_scores(
        df,
        genome_ver=genome_ver,
        pam=pam,
        scorer_config=scorer_config,
    )
    # Fill chromosome from CHOPCHOP raw genomic location when available.
    cols = {c.lower().strip(): c for c in df_raw.columns}
    genomic_col = None
    for key in ["genomic location", "genomic_location", "location", "position"]:
        if key.lower() in cols:
            genomic_col = cols[key.lower()]
            break
    if genomic_col is not None:
        extracted = (
            df_raw[genomic_col]
            .astype(str)
            .str.extract(r"(?P<chr>[^:]+):(?P<pos>\d+)")
        )
        if "chr" in extracted.columns:
            df["chr"] = extracted["chr"].fillna(default_chr if default_chr is not None else "")

    if desired_insert_pos is not None:
        df["Insert_pos"] = int(desired_insert_pos)

    return df


def _run_chopchop_from_template(loc, dist, genome_ver, pam, chopchop_config):
    if chopchop_config is None:
        chopchop_config = {}
    cmd_template = chopchop_config.get("cmd_template", "")
    if cmd_template == "":
        raise RuntimeError(
            "CHOPCHOP guide source requested, but no command template was provided. "
            "Set --chopchop_cmd_template."
        )

    chrom, pos, _strand = loc
    flank = int(chopchop_config.get("window_padding", 80))
    window_start = max(1, int(pos) - int(dist) - flank)
    window_end = int(pos) + int(dist) + flank
    target_seq = fetch_sequence_from_ensembl(
        chrom=chrom,
        start=window_start,
        end=window_end,
        strand=1,
        genome_ver=genome_ver,
        timeout=int(chopchop_config.get("ensembl_timeout", 30)),
    )
    _dump_chopchop_debug(
        chopchop_config=chopchop_config,
        source="template",
        chrom=chrom,
        pos=pos,
        window_start=window_start,
        window_end=window_end,
        target_seq=target_seq,
    )

    with tempfile.TemporaryDirectory(prefix="protospacejam_chopchop_") as temp_dir:
        fasta_path = os.path.join(temp_dir, "target.fa")
        output_file = os.path.join(temp_dir, "chopchop_output.tsv")
        with open(fasta_path, "w") as fasta_handle:
            fasta_handle.write(f">{chrom}_{pos}\n{target_seq}\n")

        format_args = {
            "fasta": fasta_path,
            "output_file": output_file,
            "output_dir": temp_dir,
            "pam": str(pam).upper(),
            "genome": genome_ver,
            "chrom": chrom,
            "pos": int(pos),
            "window_start": int(window_start),
            "window_end": int(window_end),
        }
        cmd = cmd_template.format(**format_args)
        completed = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(
                "CHOPCHOP command failed. "
                f"exit={completed.returncode}; stderr={completed.stderr.strip()}"
            )
        if not os.path.isfile(output_file):
            raise RuntimeError(
                "CHOPCHOP command completed but did not produce output file. "
                f"Expected: {output_file}"
            )

        df_raw = pd.read_csv(output_file, sep=None, engine="python")
        df = _standardize_chopchop_df(
            df_raw=df_raw,
            chr_name=chrom,
            pam=pam,
            window_start=window_start,
        )
        return _apply_external_specificity_scores(
            df,
            genome_ver=genome_ver,
            pam=pam,
            scorer_config=(chopchop_config or {}).get("scorer_config"),
        )


def _extract_job_id_from_text(text):
    if text is None:
        return None
    m = re.search(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", text, re.IGNORECASE)
    return m.group(0) if m else None


def _default_chopchop_web_payload(target_seq, genome_ver, pam, chopchop_config, gene_input=""):
    chopchop_genome = GENOME_TO_CHOPCHOP.get(genome_ver, "hg38")
    window_df = int(chopchop_config.get("window_padding", 80)) * 2 + 100
    target_type = str(chopchop_config.get("target_type", "CODING")).upper().strip()
    return {
        "opts": [
            "-J", "-BED", "-GenBank",
            "-G", chopchop_genome,
            "-filterGCmin", "10",
            "-filterGCmax", "90",
            "-t", target_type,
            "-n", "N",
            "-R", "4",
            "-P",
            "-A", "290",
            "-DF", str(window_df),
            "-a", "20",
            "-T", "1",
            "-g", "20",
            "-scoringMethod", "DOENCH_2016",
            "-f", "NN",
            "-v", "3",
            "-repairPredictions", "mESC",
            "-M", str(pam).upper(),
            "-BB", "AGGCTAGTCCGT",
        ],
        "fastaInput": "" if gene_input else target_seq,
        "geneInput": gene_input,
        "isIsoform": False,
        "forSelect": str(chopchop_config.get("for_select", "knock-in")),
    }


def _run_chopchop_web(loc, dist, genome_ver, pam, chopchop_config):
    if chopchop_config is None:
        chopchop_config = {}
    base_url = chopchop_config.get("web_base_url", "https://chopchop.cbu.uib.no").rstrip("/")
    timeout = int(chopchop_config.get("web_timeout", 180))
    poll_sec = float(chopchop_config.get("web_poll_interval", 2))

    chrom, pos, _strand = loc
    gene_input = str(chopchop_config.get("gene_input", "")).strip()
    flank = int(chopchop_config.get("window_padding", 80))
    window_start = max(1, int(pos) - int(dist) - flank)
    window_end = int(pos) + int(dist) + flank
    target_seq = ""
    if gene_input == "":
        target_seq = fetch_sequence_from_ensembl(
            chrom=chrom,
            start=window_start,
            end=window_end,
            strand=1,
            genome_ver=genome_ver,
            timeout=int(chopchop_config.get("ensembl_timeout", 30)),
        )

    payload_json = chopchop_config.get("web_payload_json", "")
    if payload_json:
        payload = json.loads(payload_json)
        payload["fastaInput"] = target_seq
        payload["geneInput"] = payload.get("geneInput", "")
        if "opts" in payload and isinstance(payload["opts"], list):
            opts = []
            i = 0
            while i < len(payload["opts"]):
                tok = payload["opts"][i]
                if tok == "-M" and i + 1 < len(payload["opts"]):
                    opts.extend(["-M", str(pam).upper()])
                    i += 2
                    continue
                if tok == "-G" and i + 1 < len(payload["opts"]):
                    opts.extend(["-G", GENOME_TO_CHOPCHOP.get(genome_ver, payload["opts"][i + 1])])
                    i += 2
                    continue
                if tok == "-t" and i + 1 < len(payload["opts"]):
                    opts.extend(["-t", str(chopchop_config.get("target_type", payload["opts"][i + 1])).upper()])
                    i += 2
                    continue
                opts.append(tok)
                i += 1
            payload["opts"] = opts
    else:
        payload = _default_chopchop_web_payload(
            target_seq=target_seq,
            genome_ver=genome_ver,
            pam=pam,
            chopchop_config=chopchop_config,
            gene_input=gene_input,
        )
    _dump_chopchop_debug(
        chopchop_config=chopchop_config,
        source="web",
        chrom=chrom,
        pos=pos,
        window_start=window_start,
        window_end=window_end,
        target_seq=target_seq,
        payload=payload,
    )

    submit_url = base_url + "/"
    req = Request(
        submit_url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "*/*",
            "Content-Type": "application/json",
            "X-Requested-With": "XMLHttpRequest",
        },
        method="POST",
    )
    with urlopen(req, timeout=timeout) as response:
        body = response.read().decode("utf-8", errors="ignore")
        header_text = str(response.headers)
        job_id = _extract_job_id_from_text(body) or _extract_job_id_from_text(header_text) or _extract_job_id_from_text(response.geturl())

    if not job_id:
        raise RuntimeError(
            "CHOPCHOP web submit succeeded but job id was not found in response. "
            "Provide --chopchop_web_payload_json copied from your browser payload if needed."
        )

    results_url = urljoin(base_url + "/", f"results/{job_id}/results.tsv")
    query_url = urljoin(base_url + "/", f"results/{job_id}/query.json")
    deadline = pd.Timestamp.utcnow().timestamp() + timeout
    last_error = None
    while pd.Timestamp.utcnow().timestamp() < deadline:
        try:
            req_results = Request(results_url, headers={"Accept": "*/*"})
            with urlopen(req_results, timeout=30) as res:
                txt = res.read().decode("utf-8", errors="ignore")
            if "Target sequence" in txt or "Rank\t" in txt:
                df_raw = pd.read_csv(StringIO(txt), sep="\t")
                try:
                    req_query = Request(query_url, headers={"Accept": "*/*"})
                    with urlopen(req_query, timeout=10):
                        pass
                except Exception:
                    pass
                win_start = window_start if gene_input == "" else 1
                df = _standardize_chopchop_df(
                    df_raw=df_raw,
                    chr_name=chrom,
                    pam=pam,
                    window_start=win_start,
                )
                return _apply_external_specificity_scores(
                    df,
                    genome_ver=genome_ver,
                    pam=pam,
                    scorer_config=(chopchop_config or {}).get("scorer_config"),
                )
        except Exception as exc:
            last_error = exc
        time.sleep(poll_sec)

    raise RuntimeError(
        f"Timed out waiting for CHOPCHOP web results ({results_url}). Last error: {last_error}"
    )


def get_chopchop_raw_results(loc, dist, genome_ver, pam, chopchop_config):
    """
    Submit a CHOPCHOP web job and return the raw results.tsv table
    with original CHOPCHOP website fields.
    """
    if chopchop_config is None:
        chopchop_config = {}
    base_url = chopchop_config.get("web_base_url", "https://chopchop.cbu.uib.no").rstrip("/")
    timeout = int(chopchop_config.get("web_timeout", 180))
    poll_sec = float(chopchop_config.get("web_poll_interval", 2))

    chrom, pos, _strand = loc
    gene_input = str(chopchop_config.get("gene_input", "")).strip()
    flank = int(chopchop_config.get("window_padding", 80))
    window_start = max(1, int(pos) - int(dist) - flank)
    window_end = int(pos) + int(dist) + flank
    target_seq = ""
    if gene_input == "":
        target_seq = fetch_sequence_from_ensembl(
            chrom=chrom,
            start=window_start,
            end=window_end,
            strand=1,
            genome_ver=genome_ver,
            timeout=int(chopchop_config.get("ensembl_timeout", 30)),
        )

    payload_json = chopchop_config.get("web_payload_json", "")
    if payload_json:
        payload = json.loads(payload_json)
        payload["fastaInput"] = target_seq
        payload["geneInput"] = payload.get("geneInput", "")
        if "opts" in payload and isinstance(payload["opts"], list):
            opts = []
            i = 0
            while i < len(payload["opts"]):
                tok = payload["opts"][i]
                if tok == "-M" and i + 1 < len(payload["opts"]):
                    opts.extend(["-M", str(pam).upper()])
                    i += 2
                    continue
                if tok == "-G" and i + 1 < len(payload["opts"]):
                    opts.extend(["-G", GENOME_TO_CHOPCHOP.get(genome_ver, payload["opts"][i + 1])])
                    i += 2
                    continue
                if tok == "-t" and i + 1 < len(payload["opts"]):
                    opts.extend(["-t", str(chopchop_config.get("target_type", payload["opts"][i + 1])).upper()])
                    i += 2
                    continue
                opts.append(tok)
                i += 1
            payload["opts"] = opts
    else:
        payload = _default_chopchop_web_payload(
            target_seq=target_seq,
            genome_ver=genome_ver,
            pam=pam,
            chopchop_config=chopchop_config,
            gene_input=gene_input,
        )
    _dump_chopchop_debug(
        chopchop_config=chopchop_config,
        source="raw_web",
        chrom=chrom,
        pos=pos,
        window_start=window_start,
        window_end=window_end,
        target_seq=target_seq,
        payload=payload,
    )

    submit_url = base_url + "/"
    req = Request(
        submit_url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "*/*",
            "Content-Type": "application/json",
            "X-Requested-With": "XMLHttpRequest",
        },
        method="POST",
    )
    with urlopen(req, timeout=timeout) as response:
        body = response.read().decode("utf-8", errors="ignore")
        header_text = str(response.headers)
        job_id = _extract_job_id_from_text(body) or _extract_job_id_from_text(header_text) or _extract_job_id_from_text(response.geturl())

    if not job_id:
        raise RuntimeError("CHOPCHOP job id not found in web submit response.")

    results_url = urljoin(base_url + "/", f"results/{job_id}/results.tsv")
    deadline = pd.Timestamp.utcnow().timestamp() + timeout
    last_error = None
    while pd.Timestamp.utcnow().timestamp() < deadline:
        try:
            req_results = Request(results_url, headers={"Accept": "*/*"})
            with urlopen(req_results, timeout=30) as res:
                txt = res.read().decode("utf-8", errors="ignore")
            if "Target sequence" in txt or "Rank\t" in txt:
                return pd.read_csv(StringIO(txt), sep="\t")
        except Exception as exc:
            last_error = exc
        time.sleep(poll_sec)

    raise RuntimeError(
        f"Timed out waiting for CHOPCHOP raw web results ({results_url}). Last error: {last_error}"
    )


def get_gRNAs_near_loc(
    loc,
    dist,
    loc2file_index,
    genome_ver,
    pam,
    guide_source="precomputed",
    chopchop_config=None,
):
    """
    input
        loc: [chr,pos,strand]
        dist: max cut to loc distance
    return:
        a dataframe of guide RNAs which cuts <[dist] to the loc, with columns in GUIDE_COLUMNS
    """
    pam = pam.upper()
    chr = loc[0]
    pos = loc[1]
    use_cut_distance_filter = False

    if str(guide_source).lower() == "chopchop":
        cmd_template = ""
        if chopchop_config is not None:
            cmd_template = chopchop_config.get("cmd_template", "")
        if cmd_template:
            df_gRNA = _run_chopchop_from_template(
                loc=loc,
                dist=dist,
                genome_ver=genome_ver,
                pam=pam,
                chopchop_config=chopchop_config,
            )
        else:
            df_gRNA = convert_chopchop_raw_to_psj(
                df_raw=get_chopchop_raw_results(
                    loc=loc,
                    dist=dist,
                    genome_ver=genome_ver,
                    pam=pam,
                    chopchop_config=chopchop_config,
                ),
                pam=pam,
                window_start=1,
                desired_insert_pos=int(pos),
                default_chr=str(chr),
                genome_ver=genome_ver,
                scorer_config=(chopchop_config or {}).get("scorer_config"),
            )
            use_cut_distance_filter = True
    else:
        if loc2file_index is None:
            return _empty_guides_df()
        chr_dict = loc2file_index[chr]
        target_files = []  # a list of file names containing gRNAs near loc
        # lookup the file
        for key in chr_dict.keys():
            file_start = key.split("-")[0]
            file_end = key.split("-")[1]
            interval = [file_start, file_end]
            if (
                in_interval(pos, interval)
                or in_interval(pos - 1000, interval)
                or in_interval(pos + 1000, interval)
            ):
                target_files.append(chr_dict[key])

        dfs = []
        for file in target_files:
            file_path = os.path.join(
                "precomputed_gRNAs",
                f"gRNAs_{pam}",
                f"gRNA_{genome_ver}",
                "gRNA.tab.gz.split.BwaMapped.scored",
                file,
            )
            df_tmp = pd.read_csv(
                file_path,
                sep="\t",
                compression="infer",
                header=None,
                names=GUIDE_COLUMNS,
            )
            dfs.append(df_tmp)
        if len(dfs) == 0:
            return _empty_guides_df()
        df_gRNA = pd.concat(dfs)

    if use_cut_distance_filter:
        tmp = df_gRNA.copy()
        tmp["_cut_pos"] = tmp.apply(lambda r: get_cut_pos(r["start"], r["end"], r["strand"]), axis=1)
        tmp["_abs_cut2ins"] = pd.to_numeric(tmp["_cut_pos"], errors="coerce").sub(float(pos)).abs()
        tmp = tmp[tmp["_abs_cut2ins"].notna()]
        tmp = tmp[tmp["_abs_cut2ins"] <= float(dist)]
        return tmp.drop(columns=["_cut_pos", "_abs_cut2ins"], errors="ignore")[GUIDE_COLUMNS]

    # subset gRNA based on strand  !ATTN: start > end when strand is '-'
    df_gRNA_on_sense = df_gRNA[(df_gRNA["strand"] == "+")]
    df_gRNA_on_antisense = df_gRNA[(df_gRNA["strand"] == "-")]

    # subset gRNAs and retain those cuts <[dist] to the loc
    df_gRNA_on_sense = df_gRNA_on_sense[
        (df_gRNA_on_sense["start"] > (pos - 17 - dist))
        & (df_gRNA_on_sense["start"] < (pos - 17 + dist))
    ]
    df_gRNA_on_antisense = df_gRNA_on_antisense[
        (df_gRNA_on_antisense["start"] > (pos + 17 - dist))
        & (df_gRNA_on_antisense["start"] < (pos + 17 + dist))
    ]

    return pd.concat([df_gRNA_on_sense, df_gRNA_on_antisense])[GUIDE_COLUMNS]


def in_interval(pos, interval):
    """
    check if pos in is interval
    input
        pos
        interval: [start,end]
    return: boolean
    """
    pos = int(pos)
    interval = list(interval)
    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    if pos >= interval[0] and pos <= interval[1]:
        # log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False


def in_interval_leftrightInclusive(pos, interval):
    """
    check if pos in is interval
    input
        pos
        interval: [start,end]
    return: boolean
    """
    pos = int(pos)
    interval = list(interval)
    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    if pos >= interval[0] and pos <= interval[1]:
        # log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False


def update_dict_count(
    key, dict
):  # update the dictionary that keeps the count of each string (key)
    if key in dict.keys():
        dict[key] += 1
    else:
        dict[key] = 1
    return dict


def get_seq(chr, start, end, strand, genome_ver):
    """
    chr
    start (1-indexed)
    end (the end position is not included
    strand 1 or -1 (str)
    """
    #print(f"fetching {chr}:{start}-{end} strand {strand}")
    chr_file_path = os.path.join("genome_files", "fa_pickle", genome_ver, f"{chr}.pk")
    log.debug(f"opening file {chr_file_path}")
    if os.path.isfile(chr_file_path):
        # read file
        chr_seqrecord = read_pickle_files(chr_file_path)
        subseq = str(chr_seqrecord.seq)[
            (start - 1) : (end - 1)
        ]  # use -1 to convert 1-index to 0-index
        if strand == "-1" or strand == -1:
            return reverse_complement(subseq)
        else:
            return subseq
    else:
        use_ensembl_fallback = os.environ.get("PROTOSPACEJAM_USE_ENSEMBL_SEQ", "1")
        if use_ensembl_fallback == "1":
            log.warning(
                f"Local genome pickle not found ({chr_file_path}); fetching region from Ensembl REST."
            )
            return fetch_sequence_from_ensembl(
                chrom=chr,
                start=start,
                end=end - 1,
                strand=strand,
                genome_ver=genome_ver,
            )
        sys.exit(f"ERROR: file not found: {chr_file_path}")


def cal_elapsed_time(starttime, endtime):
    """
    output: [elapsed_min,elapsed_sec]
    """
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    return [elapsed_min, elapsed_sec]


def read_pickle_files(file):
    if os.path.isfile(file):
        with open(file, "rb") as handle:
            mydict = pickle.load(handle)
        return mydict
    else:
        sys.exit(f"Cannot open file: {file}")


def count_ATG_at_exonEnd(ENST_info):
    """
    return a list of two items:
        count of number of ENST_IDs with ATG at the end of the exon
        list of such ENST_ID
    """
    count = 0
    list = []
    for ENST_ID in ENST_info.keys():
        my_transcript = ENST_info[ENST_ID]  # get the seq record
        if check_ATG_at_exonEnd(my_transcript):
            list.append(ENST_ID)
            count += 1
    return [count, list]


def check_ATG_at_exonEnd(my_transcript):
    """
    input: transcript object
    output: Bool
    """
    transcript_type = my_transcript.description.split("|")[1]
    if transcript_type == "protein_coding":  # only look at protein-coding transcripts
        # constructing the list of cds
        cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
        if len(cdsList) >= 1:  # has more than 1 cds
            CDS_first = cdsList[0]
            cds_len = (
                abs(CDS_first.location.start - CDS_first.location.end) + 1
            )  # for ATG to be the end of the exon, the first exon length is 3bp
            if cds_len == 3:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def get_cds_seq_in_transcript(mytranscript):
    """
    input: Bio.SeqRecord  (such as that from function fetch_ensembl_transcript)
    return the cds sequence as a string
    """
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_seqs = []
    negative_strand_flag = False
    for feat in mytranscript.features:
        if feat.type == "cds":
            st = feat.location.start.position
            en = feat.location.end.position
            if feat.strand == -1:  # neg strand
                cds_seqs.append(
                    str(Seq(wholeSeq_rc[st:en]).reverse_complement())
                )  # the coord are respective to the revcom of the retrieved seq (weird)
                negative_strand_flag = True
            else:  # pos strand
                cds_seqs.append(wholeSeq[st:en])
    if negative_strand_flag == True:
        return "".join(cds_seqs[::-1])
    else:
        return "".join(cds_seqs)


def get_cds_seqNflank(transcriptObj, which_cds, cds_flank_len):
    """
    returns the n-th cds and its flanking sequences
    the returned sequence will be in the coding strand
    """
    total_cds_num = transcriptObj.num_cds
    if which_cds > total_cds_num:
        raise ValueError(
            f"Requested coding exon is {which_cds}, but the transcript only has {total_cds_num} coding exons"
        )

    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in transcriptObj.features if feat.type == "cds"]
    strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    which_cds -= 1

    if strand == 1:
        which_cds = which_cds
    else:
        which_cds = len(cds_list) - which_cds - 1

    # get cds seq with flank
    target_cds = cds_list[which_cds]
    cds_st = target_cds.location.start.position
    cds_en = target_cds.location.end.position
    if strand == 1:
        cds_seq = wholeSeq[cds_st:cds_en]
        Lflank_st = max([i for i in [cds_st - cds_flank_len, 0] if i >= 0])
        RFlank_en = min(
            [i for i in [cds_en + cds_flank_len, len(wholeSeq)] if i <= len(wholeSeq)]
        )
        Lflank = wholeSeq[Lflank_st:cds_st]
        Rflank = wholeSeq[cds_en:RFlank_en]
    else:
        cds_seq = str(
            Seq(wholeSeq_rc[cds_st:cds_en]).reverse_complement()
        )  # the coord are respective to the revcom of the retrieved seq (weird)
        Lflank_st = max([i for i in [cds_st - cds_flank_len, 0] if i >= 0])
        RFlank_en = min(
            [
                i
                for i in [cds_en + cds_flank_len, len(wholeSeq_rc)]
                if i <= len(wholeSeq_rc)
            ]
        )
        Lflank = str(Seq(wholeSeq_rc[Lflank_st:cds_st]).reverse_complement())
        Rflank = str(Seq(wholeSeq_rc[cds_en:RFlank_en]).reverse_complement())
        Lflank, Rflank = Rflank, Lflank  # for -1 strand, the Lflank is the Rflank

    return {"cds_seq": cds_seq, "Lflank": Lflank, "Rflank": Rflank}


def get_exon_concat_with_flank(transcriptObj, which_cds, cds_flank_len):
    """
    get the nth cds segment from the transcript and concat with flanks
    the returned sequence will be in the coding strand
    """
    cds_seq_N_flank = get_cds_seqNflank(
        transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len
    )
    cds_with_flank = "".join(
        [
            cds_seq_N_flank["Lflank"],
            cds_seq_N_flank["cds_seq"],
            cds_seq_N_flank["Rflank"],
        ]
    )
    return cds_with_flank


def get_cutsite_in_gene(transcriptObj, listOfgRNAObj, which_cds, cds_flank_len):
    """
    find the coordinate (respective to the gene) of the cutsite in each gRNA

    Parameters
    ----------
    transcriptObj
    listOfgRNAObj
    which_cds: the i-th cds used to extract the gRNA
    cds_flank_len: the len of flank added to the cds when extracting gRNA sequences

    Returns
    -------
    list of gRNA object (updated with the cutsite)
    """
    cds_w_flank = get_exon_concat_with_flank(
        transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len
    )

    # find cds start and end
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_list = [feat for feat in transcriptObj.features if feat.type == "cds"]
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    cds_human_readable = which_cds
    which_cds -= 1
    if cds_strand == 1:
        which_cds = which_cds
    else:
        which_cds = len(cds_list) - which_cds - 1

    target_cds = cds_list[which_cds]
    cds_st = target_cds.location.start.position
    cds_en = target_cds.location.end.position
    cds_len = cds_en - cds_st

    # update the cutsite for each gRNA
    for idx, gRNAObj in enumerate(listOfgRNAObj):
        # print(f"{gRNAObj.protospacer} {gRNAObj.pam} {gRNAObj.g_strand} {gRNAObj.g_st} {gRNAObj.g_en}")

        # find gRNA cut site in cds
        gRNACut_in_cds = (
            gRNAObj.g_st - cds_flank_len + 17
        )  # the num of bp before cutsite

        if gRNAObj.g_strand == "-":
            gRNACut_in_cds = (len(cds_w_flank) - 2 * cds_flank_len) - gRNACut_in_cds

        # find gRNA cut site in gene, this is relative to the reference_left_index
        gRNACut_in_gene = gRNACut_in_cds + cds_st
        if (
            cds_strand == -1
        ):  # RNACut_in_cds needs to be from the right if strand == -1 (making it relative to the reference_left_index)
            tmp_gRNACut_in_cds = (len(cds_w_flank) - 2 * cds_flank_len) - gRNACut_in_cds
            gRNACut_in_gene = tmp_gRNACut_in_cds + cds_st

        # find gRNA cut site in chromosome
        gRNACut_in_chr = (
            gRNACut_in_gene + transcriptObj.annotations["reference_left_index"]
        )

        listOfgRNAObj[idx].gRNACut_in_cds = gRNACut_in_cds
        listOfgRNAObj[idx].gRNACut_in_gene = gRNACut_in_gene
        listOfgRNAObj[idx].gRNACut_in_chr = gRNACut_in_chr
        listOfgRNAObj[idx].cds_len = cds_len
        listOfgRNAObj[idx].cds = cds_human_readable
        listOfgRNAObj[idx].Ensemble_ID = transcriptObj.id
        listOfgRNAObj[idx].Ensemble_ref = transcriptObj.annotations["reference_species"]
        listOfgRNAObj[idx].Ensemble_chr = transcriptObj.annotations[
            "reference_chromosome_number"
        ]
        listOfgRNAObj[idx].Ensemble_chr_left_idx = transcriptObj.annotations[
            "reference_left_index"
        ]
        listOfgRNAObj[idx].Ensemble_chr_right_idx = transcriptObj.annotations[
            "reference_right_index"
        ]
        listOfgRNAObj[idx].Ensemble_transcript_strand = transcriptObj.annotations[
            "transcript_strand"
        ]

    return listOfgRNAObj


def get_HDR_flank(transcriptObj, listOfgRNAObj, HDR_flank_len):
    """
    Parameters
    ----------
    transcriptObj
    listOfgRNAObj
    HDR_flank_len: [int] the desired length of the HDR length
    which_cds: the i-th cds used to extract the gRNA (this is used for adjusting the frame)
    adjust_frame: [bool] True/False

    Returns
    -------
    list of gRNA object (updated with the HDR flanks, and whether if the HDR flank is shorter than expected, flags: HDR_Rflank_short and HDR_Lflank_short)
    """
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        listOfgRNAObj[idx].HDR_Lflank_short = 0
        listOfgRNAObj[idx].HDR_Rflank_short = 0
        gRNACut_in_gene = int(gRNAObj.gRNACut_in_gene)
        HDR_Lflank = ""
        HDR_Rflank = ""
        if cds_strand == 1:
            L_st = gRNACut_in_gene - HDR_flank_len
            L_en = gRNACut_in_gene
            R_st = gRNACut_in_gene
            R_en = gRNACut_in_gene + HDR_flank_len

            if L_st < 0:  # check out of bounds
                L_st = 0
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_st < 0:
                R_st = 0
                listOfgRNAObj[idx].HDR_Rflank_short = 1
            if L_en >= gene_len:
                L_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_en >= gene_len:
                R_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Rflank_short = 1

            HDR_Lflank = wholeSeq[L_st:L_en]
            HDR_Rflank = wholeSeq[R_st:R_en]

        else:  # strand = -1
            L_st = gRNACut_in_gene
            L_en = gRNACut_in_gene + HDR_flank_len
            R_st = gRNACut_in_gene - HDR_flank_len
            R_en = gRNACut_in_gene

            if L_st < 0:  # check out of bounds
                L_st = 0
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_st < 0:
                R_st = 0
                listOfgRNAObj[idx].HDR_Rflank_short = 1
            if L_en >= gene_len:
                L_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_en >= gene_len:
                R_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Rflank_short = 1

            HDR_Lflank = str(Seq(wholeSeq_rc[L_st:L_en]).reverse_complement())
            HDR_Rflank = str(
                Seq(wholeSeq_rc[R_st:R_en]).reverse_complement()
            )  # the coord are respective to the revcom of the retrieved seq (weird)

        listOfgRNAObj[idx].HDR_Lflank = HDR_Lflank
        listOfgRNAObj[idx].HDR_Rflank = HDR_Rflank

    return listOfgRNAObj


def calculate_offset_cutsite(
    gRNACut_in_cds, leading_nonTriplet
):  # calculate the offset so the edit will be in frame
    """
    Parameters
    ----------
    leading_nonTriplet

    Returns
    -------
    adjusted L_en
    """
    adj = 0
    if gRNACut_in_cds == 0:
        if leading_nonTriplet == 0:
            adj = 0
        elif leading_nonTriplet == 1:
            adj = 1
        elif leading_nonTriplet == 2:
            adj = -1
    elif gRNACut_in_cds == 1:
        if leading_nonTriplet == 0:
            adj = -1
        elif leading_nonTriplet == 1:
            adj = 0
        elif leading_nonTriplet == 2:
            adj = 1
    elif gRNACut_in_cds == 2:
        if leading_nonTriplet == 0:
            adj = 1
        elif leading_nonTriplet == 1:
            adj = -1
        elif leading_nonTriplet == 2:
            adj = 0
    elif gRNACut_in_cds >= 3:
        trailing_nonTriplet = (gRNACut_in_cds - leading_nonTriplet) % 3
        if trailing_nonTriplet == 0:
            adj = 0
        elif trailing_nonTriplet == 1:
            adj = -1
        elif trailing_nonTriplet == 2:
            adj = 1

    return adj


def nudge_cutsite_inframe(transcriptObj, listOfgRNAObj, HDR_flank_len, which_cds):
    """
    nudge the cutsite so that it is in frame (although it doesn't reflect the true cutsite now)
    Parameters
    ----------
    L_en: Lflank end
    leading_nonTriplet

    Returns
    -------
    adjusted L_en
    """
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    # calculate the number of leading nonTriplet at the 5' of the cds/exon
    frame = get_cds_frame(transcriptObj, which_cds)
    if frame == 1:
        leading_nonTriplet = 0
    if frame == 2:
        leading_nonTriplet = 1
    if frame == 3:
        leading_nonTriplet = 2

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        gRNACut_in_cds = int(gRNAObj.gRNACut_in_cds)
        adj = calculate_offset_cutsite(gRNACut_in_cds, leading_nonTriplet)
        # nudge cutsite (respective to the gene)
        if cds_strand == -1:
            adj = 0 - adj
        listOfgRNAObj[idx].gRNACut_in_gene = listOfgRNAObj[idx].gRNACut_in_gene + adj
        listOfgRNAObj[idx].cds_leading_nonTriplet = leading_nonTriplet
    return listOfgRNAObj


def get_cds_frame(mytranscript, which_cds):
    """
    for the n-th cds, calculate the frame (using coding sequence from 1 to n-1 th cds)
    """
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in mytranscript.features if feat.type == "cds"]
    strand = list(
        set([feat.strand for feat in mytranscript.features if feat.type == "cds"])
    )[0]
    which_cds -= 1
    cds_until = ""

    if strand == 1:
        which_cds = which_cds
        for i in range(0, which_cds):
            st = cds_list[i].location.start.position
            en = cds_list[i].location.end.position
            cds_until = cds_until + wholeSeq[st:en]
    else:
        which_cds = len(cds_list) - which_cds - 1
        for i in range(len(cds_list) - 1, which_cds, -1):
            st = cds_list[i].location.start.position
            en = cds_list[i].location.end.position
            cds_until = cds_until + str(Seq(wholeSeq_rc[st:en]).reverse_complement())

    # calculate frame
    frame = len(cds_until) % 3
    if frame == 0:
        frame = 1
    elif frame == 1:
        frame = 3
    elif frame == 2:
        frame = 2
    return frame

if __name__ == "__main__":
    import doctest

    doctest.testmod()
