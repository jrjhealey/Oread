#!/usr/bin/env pythonw

# TODO:
# Consider incorporating a contig reorder?

import os
import sys
import logging
import tempfile

from gooey import Gooey, GooeyParser

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline


logging.basicConfig(format="[%(asctime)s] %(levelname)-8s \n -> %(message)s",
                    level=logging.NOTSET, datefmt='%m/%d/%Y %I:%M:%S%p')
logger = logging.getLogger(__name__)


@Gooey(program_name="Oread - the companion to Artemis.",
       error_color="#cb4b16",
       header_bg_color="#2aa198",
       default_size=(750, 550))
def get_args():
    """Parse command line arguments"""
    desc = (" /ˈɔːrɪad/\n"
            "Create genome comparison files for use with Artemis/ACT.")


    try:
        parser = GooeyParser(description=desc, prog="oread")
        group = parser.add_argument_group()
        group.add_argument(
            "-s",
            "--subject",
            action="store",
            metavar="Subject",
            required=True,
            widget="FileChooser",
            help="Subject sequence file (.fasta)."
        )
        group.add_argument(
            "-q",
            "--query",
            action="store",
            metavar="Query",
            required=True,
            widget="FileChooser",
            help="Query sequence file (.fasta)."
        )
        group.add_argument(
            "-o",
            "--outdir",
            action="store",
            metavar="Outdir comparison folder",
            widget="DirChooser",
            help="Your outdir folder."
        )
        # It seems running blast in this mode doesnt allow for multiple threads...
        # group.add_argument(
        #     "-c",
        #     "--cpus",
        #     metavar="Number of threads for BLAST",
        #     widget="Dropdown",
        #     default=str(round(os.cpu_count()/2)),
        #     choices=[str(x) for x in range(1, os.cpu_count()+1)],
        #     help="Number of threads for BLAST to use."
        # )
        group.add_argument(
            "-k",
            "--keep_temp",
            metavar="Keep any temporary files produced.",
            widget="Dropdown",
            choices=["True", "False"],
            default="False",
            help="Retain temporary files produced during the run."
        )
        group.add_argument(
            "-t",
            "--task",
            metavar="Type of BLAST search.",
            widget="Dropdown",
            choices=["megablast", "dc-megablast", "blastn"],
            default="blastn",
            help="What type of BLAST to run:\n"
                 " - Megablast for very similar sequences (e.g. sequencing errors)\n"
                 " - dc-megablast typically for inter-species comparison\n"
                 " - 'traditional' blastn for more dissimilar sequences."
        )

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def format_outfile(args):
    """Synthesise an outdir filename for the comparison file"""
    return os.path.join(args.outdir, "{}_vs_{}.act".format(
            os.path.basename(os.path.splitext(args.subject)[0]),
            os.path.basename(os.path.splitext(args.query)[0])))


def basename(string):
    """Abstraction/wrapper for file basenames from os.path"""
    return os.path.splitext(os.path.abspath(string))[0]


def main():
    """Main method for creation of Artemis comparison files"""
    logger.info("Launching {}...".format(__file__))

    args = get_args()

# Synthesise the outdir name, if the directory is specified, else synthesise both

    if not args.outdir:
        args.outdir = os.path.dirname(os.path.abspath(args.subject))
        outfile = format_outfile(args)
    if args.outdir:
        outfile = format_outfile(args)

        logger.info("No output directory specified, defaulting to: {}".format(args.outdir))
        logger.info("Comparison file will be stored at {}".format(outfile))

# ACT cannot handle multiblast files (i.e. files generated from multifastas etc.
# Therefore, a temporary concatenated file is needed.
# Check if this is needed, and make temporary files as required.

    intermediate_query_needed = False
    intermediate_subj_needed = False
    try:
        SeqIO.read(args.query, "fasta")
    except ValueError:
        logger.warning("More than one sequence detected in the query file.")
        intermediate_query_needed = True
    try:
        SeqIO.read(args.subject, "fasta")
    except ValueError:
        logger.warning("More than one sequence detected in the subject file.")
        intermediate_subj_needed = True
    finally:
        logger.warning("ACT Comparison files require contiguous sequences. "
                       "Temporary files will be created, but the comparison file "
                       "will still work with the unconcatenated sequences inside ACT.")

    if intermediate_query_needed or intermediate_subj_needed:
        tempfile.tempdir = args.outdir
        if intermediate_query_needed:
            logger.info("Creating intermediate files...")
            logger.info("Storing temporary files in {}".format(tempfile.tempdir))

            with tempfile.NamedTemporaryFile(mode="w",prefix=basename(args.query),
                                             suffix=".fa", delete=False) as iqry:
                logger.info("Temporary query file is {}".format(iqry.name))

                iqry.write(">{}_intermediate_concatenation\n".format(basename(args.query)))
                for record in SeqIO.parse(args.query, "fasta"):
                    iqry.write(str(record.seq.rstrip("\n")))

                iqry.seek(0)
                    # iqry.write("\n".join(wrap(data['Sequence'], 60)) + "\n")
                    # snippet for wrapping if desired later.
                args.query = iqry.name

        if intermediate_subj_needed:
            with tempfile.NamedTemporaryFile(mode="w", prefix=basename(args.subject),
                                             suffix=".fa", delete=False) as isub:
                logger.info("Temporary subject file is {}".format(isub.name))

                isub.write(">{}_intermediate_concatenation\n".format(basename(args.subject)))
                for record in SeqIO.parse(args.subject, "fasta"):
                    isub.write(str(record.seq.rstrip("\n")))

                isub.seek(0)

                args.subject = isub.name


    blast = NcbiblastnCommandline(
        cmd="blastn",
        task=args.task,
        query=args.query,
        subject=args.subject,
        outfmt=6,
        out=outfile
    )

    logger.info("Running BLASTn, as follows:")
    logger.info(str(blast))
    stdout, stderr = blast()
    if stderr:
        sys.stderr.write(stderr)

    if intermediate_query_needed:
        os.unlink(args.query)
    if intermediate_subj_needed:
        os.unlink(args.query)

if __name__ == "__main__":
    main()
