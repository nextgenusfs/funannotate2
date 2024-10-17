#!/usr/bin/env python3

import sys
import os
import argparse
from .__version__ import __version__
from .help_formatter import MyParser, MyHelpFormatter

from .annotate import annotate

# from .update import update
from .train import train
from .predict import predict
from .clean import clean
from .install import install


def main():
    args = parse_args(sys.argv[1:])
    if args.subparser_name == "install":
        install(args)
    elif args.subparser_name == "clean":
        clean(args)
    elif args.subparser_name == "predict":
        predict(args)
    elif args.subparser_name == "train":
        train(args)
    elif args.subparser_name == "annotate":
        annotate(args)
    # elif args.subparser_name == 'update':
    #    update(args)


def parse_args(args):
    description = "Funannotate: eukaryotic genome annotation pipeline"
    parser = MyParser(
        description=description, formatter_class=MyHelpFormatter, add_help=False
    )
    subparsers = parser.add_subparsers(title="Commands", dest="subparser_name")
    # add subtools here
    install_subparser(subparsers)
    clean_subparser(subparsers)
    train_subparser(subparsers)
    predict_subparser(subparsers)
    annotate_subparser(subparsers)
    # update_subparser(subparsers)

    help_args = parser.add_argument_group("Help")
    help_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )
    help_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(args)


def install_subparser(subparsers):
    group = subparsers.add_parser(
        "install",
        description="Install funannotate2 database files, $FUNANNOTATE2_DB env variable is required to be set.|n",
        help="Install funannotate2 database files.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )
    required_args = group.add_argument_group("Required arguments")
    required_args.add_argument(
        "-d",
        "--db",
        required=True,
        nargs="+",
        choices=[
            "all",
            "merops",
            "uniprot",
            "dbCAN",
            "pfam",
            "go",
            "mibig",
            "interpro",
            "gene2product",
        ],
        help="Databases to install [all,merops,uniprot,dbCAN,pfam,go,mibig,interpro,gene2product]",
        metavar="",
    )
    optional_args = group.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-w",
        "--wget",
        action="store_true",
        help="Use wget for downloading",
    )
    optional_args.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force re-download/re-install of all databases",
    )
    optional_args.add_argument(
        "-u",
        "--update",
        action="store_true",
        help="Update databases if change detected",
    )
    other_args = group.add_argument_group("Other arguments")
    other_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    other_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )


def predict_subparser(subparsers):
    group = subparsers.add_parser(
        "predict",
        description="""Gene model prediction from automated training of ab initio gene prediction algorithms:|n\
""",
        help="Predict primary gene models in a eukaryotic genome.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )
    required_args = group.add_argument_group("Required arguments")
    required_args.add_argument(
        "-i",
        "--input-dir",
        required=False,
        dest="input_dir",
        help="funannotate2 output directory",
        metavar="",
    )
    required_args.add_argument(
        "-f",
        "--fasta",
        required=False,
        help="genome in FASTA format (softmasked repeats)",
        metavar="",
    )
    required_args.add_argument(
        "-o", "--out", required=False, help="Output folder name", metavar=""
    )
    required_args.add_argument(
        "-p",
        "--params",
        "--pretrained",
        dest="params",
        required=False,
        help="Params.json or pretrained species slug. `funannotate2 info` to see pretrained species",
        metavar="",
    )
    required_args.add_argument(
        "-s",
        "--species",
        required=False,
        help='Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"',
        metavar="",
    )
    optional_args = group.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-st",
        "--strain",
        required=False,
        help="Strain/isolate name (e.g. Af293)",
        metavar="",
    )
    optional_args.add_argument(
        "-w",
        "--weights",
        required=False,
        nargs="+",
        help="Gene predictors and weights",
        metavar="",
    )
    optional_args.add_argument(
        "-m",
        "--consensus-method",
        required=False,
        choices=["evm", "gfftk"],
        dest="consensus",
        help="Consensus model generation method [evm,gfftk]",
        default="gfftk",
        metavar="",
    )
    optional_args.add_argument(
        "-ps",
        "--proteins",
        required=False,
        help="Proteins to use for evidence",
        nargs="+",
        metavar="",
    )
    optional_args.add_argument(
        "-ts",
        "--transcripts",
        required=False,
        help="Transcripts to use for evidence",
        nargs="+",
        metavar="",
    )
    optional_args.add_argument(
        "-c",
        "--cpus",
        help="Number of cpus/threads to use",
        type=int,
        default=2,
        metavar="",
    )
    optional_args.add_argument(
        "-mi",
        "--max-intron",
        dest="max_intron",
        help="Maximum intron length",
        type=int,
        default=3000,
        metavar="",
    )
    optional_args.add_argument(
        "-hl",
        "--header-len",
        default=100,
        dest="header_length",
        type=int,
        help="Max length for fasta headers",
        metavar="",
    )
    optional_args.add_argument(
        "-l",
        "--locus-tag",
        dest="name",
        default="FUN_",
        help="Locus tag for genes, perhaps assigned by NCBI, eg. VC83",
        metavar="",
    )
    optional_args.add_argument(
        "-n",
        "--numbering",
        default=1,
        type=int,
        help="Specify start of gene numbering",
        metavar="",
    )
    optional_args.add_argument(
        "-pt",
        "--pretrained-species",
        dest="pretrained_species",
        required=False,
        help="Use pretrained parameters for prediction, ie values in funannotate2 species",
        metavar="",
    )
    optional_args.add_argument(
        "--tmpdir", default="/tmp", help="volume to write tmp files", metavar=""
    )
    other_args = group.add_argument_group("Other arguments")
    other_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    other_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )


def train_subparser(subparsers):
    group = subparsers.add_parser(
        "train",
        description="Gene model training from automated training of ab initio gene prediction algorithms:|n",
        help="Train ab intio gene prediction algorithms.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )
    required_args = group.add_argument_group("Required arguments")
    required_args.add_argument(
        "-f", "--fasta", required=True, help="genome in FASTA format", metavar=""
    )
    required_args.add_argument(
        "-s",
        "--species",
        required=True,
        help='Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"',
        metavar="",
    )
    required_args.add_argument(
        "-o", "--out", required=True, help="Output folder name", metavar=""
    )
    optional_args = group.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-t",
        "--training-set",
        dest="training_set",
        help="Training set to use in GFF3 format",
        metavar="",
    )
    optional_args.add_argument("--strain", help="Strain/isolate name", metavar="")
    optional_args.add_argument(
        "--cpus", default=2, type=int, help="Number of CPUs to use", metavar=""
    )
    optional_args.add_argument(
        "--optimize-augustus",
        action="store_true",
        dest="optimize_augustus",
        help="Run Augustus mediated optimized training (not recommended)",
    )
    optional_args.add_argument(
        "--header-len",
        default=100,
        dest="header_length",
        type=int,
        help="Max length for fasta headers",
        metavar="",
    )
    other_args = group.add_argument_group("Other arguments")
    other_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    other_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )


def clean_subparser(subparsers):
    group = subparsers.add_parser(
        "clean",
        description="The script sorts contigs by size, starting with shortest contigs it uses minimap2 to find contigs duplicated elsewhere, and then removes duplicated contigs.",
        help="Find and remove duplicated contigs, sort by size, rename headers.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )

    required_args = group.add_argument_group("Required arguments")
    required_args.add_argument(
        "-f", "--fasta", required=True, help="genome in FASTA format", metavar=""
    )
    required_args.add_argument(
        "-o", "--out", required=True, help="Cleaned genome output (FASTA)", metavar=""
    )

    optional_args = group.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-p",
        "--pident",
        type=int,
        default=95,
        help="percent identity of contig",
        metavar="",
    )
    optional_args.add_argument(
        "-c", "--cov", type=int, default=95, help="coverage of contig", metavar=""
    )
    optional_args.add_argument(
        "-m",
        "--minlen",
        type=int,
        default=500,
        help="Minimum length of contig",
        metavar="",
    )
    optional_args.add_argument(
        "-r",
        "--rename",
        help="Rename contigs largest to smallest with this basename, ie scaffold_",
        metavar="",
    )
    optional_args.add_argument(
        "--cpus", default=2, type=int, help="Number of CPUs to use", metavar=""
    )
    optional_args.add_argument("--tmpdir", help="TMP directory to use", metavar="")
    optional_args.add_argument(
        "--exhaustive",
        action="store_true",
        help="Compute every contig, else stop at N50",
    )
    optional_args.add_argument(
        "--logfile", default=False, help="Write logs to file", metavar=""
    )
    optional_args.add_argument("--debug", action="store_true", help="Debug the output")

    other_args = group.add_argument_group("Other arguments")
    other_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    other_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )


def annotate_subparser(subparsers):
    group = subparsers.add_parser(
        "annotate",
        description="Add functional annotate to gene models|n",
        help="Add functional annotation to gene models.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )
    required_args = group.add_argument_group("Required arguments")
    required_args.add_argument(
        "-i",
        "--input-dir",
        required=False,
        dest="input_dir",
        help="funannotate2 output directory",
        metavar="",
    )
    required_args.add_argument(
        "-f", "--fasta", required=False, help="genome in FASTA format", metavar=""
    )
    required_args.add_argument(
        "-t",
        "--tbl",
        required=False,
        help="genome annotation in TBL format",
        metavar="",
    )
    required_args.add_argument(
        "-g",
        "--gff3",
        required=False,
        help="genome annotation in GFF3 format",
        metavar="",
    )
    required_args.add_argument(
        "-o", "--out", required=False, help="Output folder name", metavar=""
    )
    optional_args = group.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-s",
        "--species",
        help='Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"',
        metavar="",
    )
    optional_args.add_argument("--strain", help="Strain/isolate name", metavar="")
    optional_args.add_argument(
        "--cpus", default=2, type=int, help="Number of CPUs to use", metavar=""
    )
    optional_args.add_argument(
        "--tmpdir", default="/tmp", help="volume to write tmp files", metavar=""
    )
    other_args = group.add_argument_group("Other arguments")
    other_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    other_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )


if __name__ == "__main__":
    main()
