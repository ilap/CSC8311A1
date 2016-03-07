import sys
import tempfile
from Bio import SearchIO

from subprocess import Popen
from subprocess import PIPE


class GuideRNABioUtils:

    def __init__(self):

        self._arch = sys.platform

        if self._arch.startswith("linux"):
            # Blat executables only run on x64 Linux
            if sys.maxsize > 2**32:
                self._blat = "blat_linux_x64"
            else:
                raise Exception('BLAT only runs on 64bit Linux')
        elif self._arch == "darwin":
            self._blat = "blat_macos"
            # Mac OS X
        else:
            raise Exception('Platform {} is not currently supported. '
                            'Exiting.'.format(self._arch))

        # Initial BLAT parameters

    # Run Blat on Species model
    def run_blat(self, species_file, target_file):

        print "RUN_BLAT..."
        # print blast_parser

        temp_file = tempfile.NamedTemporaryFile()
        output = temp_file.name

        out_format = "psl"

        blat_cmd = "./utils/{0} -out={1} -t=DNA -q=DNA {2} {3} {4}".format(
            self._blat, out_format, species_file, target_file, output)

        print "BLAT CMD", blat_cmd
        blat = Popen(blat_cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=(
            sys.platform != "win32"))

        stdout, stderr = blat.communicate()
        print "STDOUT", stdout
        print "STDERR", stderr
        if stderr:
            raise Exception("blat had error during run.")

        blat_qresult = SearchIO.read(output, 'blat-{}'.format(out_format))
        print blat_qresult

        # http://biopython.org/DIST/docs/tutorial/Tutorial.html
        for hit in blat_qresult:

            print "HIT", dir(hit), hit.seq_len
            for hsp in hit:
                print "HSP", dir(hsp)
                print "....", hsp.hit_start, hsp.hit_end, hsp.hit_strand, \
                    hsp.query_strand, hsp.query_is_protein
                print hsp.hit_end - hsp.hit_start
