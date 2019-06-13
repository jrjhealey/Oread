#!/usr/bin/env pythonw

import os
import sys
import logging
import subprocess

from argparse import ArgumentParser
from gooey import Gooey, GooeyParser

from Bio.Blast.Applications import NcbiblastxCommandline

logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger(__name__)


@Gooey(program_name="Artemisia")
def get_args():
    """Parse command line arguments"""
    desc = """Create genome comparison files for use with Artemis/ACT."""

    try:
        parser = GooeyParser(description=desc, prog="artemisia")
        parser.add_argument(
            "-s",
            "--subject",
            action="store",
            metavar="Subject",
            required=True,
            widget="FileChooser",
            help="Subject sequence file (.fasta)."
        )
        parser.add_argument(
            "-q",
            "--query",
            action="store",
            metavar="Query",
            required=True,
            widget="FileChooser",
            help="Query sequence file (.fasta)."
        )
        parser.add_argument(
            "-o",
            "--output",
            action="store",
            metavar="Output comparison file",
            widget="FileChooser",
            help="Your output comparison filename (BLAST tabular)."
        )

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    """Main method for creation of Artemis comparison files"""

    args = get_args()

    if not args.output:
        args.output = os.path.join(os.path.dirname(os.path.abspath(args.subject)),
                                   "{}_vs_{}.act".format(
                                       os.path.basename(os.path.splitext(args.subject)[0]),
                                       os.path.basename(os.path.splitext(args.query)[0])))

        logger.info("No output file specified, defaulting to: " + args.output)

    blast = NcbiblastxCommandline(
        cmd="blastn",
        query=args.query,
        subject=args.subject,
        outfmt=6,
        out=args.output
    )

    logger.info("Running BLASTn, as follows:")
    logger.info(str(blast))
    stdout, stderr = blast()




if __name__ == "__main__":
    main()
