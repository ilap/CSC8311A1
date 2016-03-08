import sys
import tempfile
from Bio import SearchIO

from subprocess import Popen
from subprocess import PIPE

from .models import *


class GuideRNAManager:

    # Global private attributes
    _initialised = False

    # Instances of the selected model objects
    _species = None
    _target = None
    _pam = None

    # New models for the objects created by BLAT and search gRNA.
    _target_hit_model = None
    _guide_rna_model = None

    # List of the target gene/sequence found in the genome.
    _hits = None

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
    def _run_blat(self):

        if not self._initialised:
            return None

        species_file = self._species.fasta_file.path
        target_file = self._target.sequence_file.path

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
        # print "STDOUT", stdout
        # print "STDERR", stderr
        if stderr:
            raise Exception("blat had error during run.")

        blat_qresult = SearchIO.read(output, 'blat-{}'.format(out_format))
        print blat_qresult

        # http://biopython.org/DIST/docs/tutorial/Tutorial.html
        # Create TargetHit models from the BLAT results.
        for hit in blat_qresult:

            # BLAT hit length means, number of HSPs in Hit, if it's biffer
            # than 1 then raise and exception.
            # TODO: if there is more than one HSP then raise an error and
            # handle it properly
            if len(hit) > 1:
                raise Exception("There is more than one HSP in BLAT result.")

            # For debug purpose
            print "HIT ID", hit.id
            print "HIT SEQ LEN", hit.seq_len

            '''print "HIT LEN means HSP", len(hit)
            for hsp in hit:
                print "HIT Strand", hsp[0].hit_strand
                print "HIT Query Strand", hsp[0].query_strand
                print "HSP hit_end", hsp.hit_end
                print "HSP hit_gap_num", hsp.hit_gap_num
                print "HSP hit_gapopen_num", hsp.hit_gapopen_num
                print "HSP hit_span_all", hsp.hit_span_all
                print "HSP hit_start", hsp.hit_start
                print "HSP hit_start_all", hsp.hit_start_all
                print "HSP match_num", hsp.match_num
                print "HSP mismatch_num", hsp.mismatch_num
                print "HSP match_rep_num", hsp.match_rep_num
                print "HSP n_num", hsp.n_num
                print "HSP query_end", hsp.query_end
                print "HSP query_gap_num", hsp.query_gap_num
                print "HSP query_gapopen_num", hsp.query_gapopen_num
                print "HSP query_span_all", hsp.query_span_all
                print "HSP query_start", hsp.query_start
                print "HSP query_start_all", hsp.query_start_all
                print "LEN:", hsp.hit_end - hsp.hit_start
                print "HSP SCORE:", hsp.score'''

            length = hit[0].hit_end - hit[0].hit_start
            print "LENGTH", float(length)
            # Strand on HSPFragment, meaning: 0: 5-3, 1: 3-5
            strand = hit[0][0].hit_strand

            # hsp = hit[0]
            # hsp.match_num/length
            score = hit[0].match_num/float(length)
            print "SCORE", score

            hit_model = TargetHit.objects.create(species=self._species,
                                                 target=self._target,
                                                 position=hit[0].hit_start,
                                                 strand=strand,
                                                 score=score,
                                                 length=length)

            hit_model.save()

            return True

    def initialise_run(self, request):

        if self._initialised:
            raise Exception("The GuideRNABioutisl should be singleton")
        species_id = request.POST['species']
        target = request.POST['target']
        pam_id = request.POST['pam']
        upstream = request.POST['upstream']
        downstream = request.POST['downstream']
        is_nickase = request.POST['is_nickase']

        # Steps are the following:
        # 1. Get the Species object
        # TODO: Restructure PAM based on Nucleases are used
        # 2. get the PAM Object or all if it's -1
        # TODO: Make available for multiple genome and/or targets
        #  3. Create Target object

        # 1.
        self._species = Species.objects.get(pk=species_id)

        # 2.
        # -1 means all PAM should be used for search gRNA

        if (pam_id<0):
            self._pam = PAM.objects.all()
        else:
            self._pam = PAM.objects.get(pk=pam_id)

        # 3. Create target model from the seqeuence
        self._target = Target()
        self._target.init_target(target, upstream, downstream)

        self._initialised = True
        return self

    # TODO: Only searches for guide RNA in the target sequence +- offset.
    def search_grna(self):
        if not self._initialised:
            return None

        # After the GuideRNAManager object is initialised
        # 1. Run BLAT on # target against species genome
        # 2. Run search grna on target hits found in genome.
        # 3. Visualise the result or run

        blat_results = self._run_blat()
        if not blat_results:
            raise Exception("Run BLAT is failed.")

        query_file = self._species.fasta_file.path
        target_file = self._target.sequence_file.path
        print "QUERY", query_file
        print "TARGET", target_file

        # TODO: Only one FASTA record is allowed from the query and target.
        query_record = SeqIO.read(query_file, "fasta")
        target_record = SeqIO.read(target_file, "fasta")

        print query_record.id, target_record.id
        print repr(query_record.seq)
        print len(query_record)

        print target_record.id
        print repr(target_record.seq)
        print len(target_record)

        '''output_handle = open("/tmp/example.fasta", "w")
        output_handle.write(str(query_record.seq[:]))
        output_handle.close()'''

        # TODO: Only searches in the target hits of the genome
        for hit in TargetHit.objects.filter(species=self._species,
                                             target=self._target):
            print "HIT STRAND", hit.strand
            target_start = hit.position
            target_end = hit.position+hit.length
            query_sequence = query_record.seq[target_start:target_end]

            print "HINT", hit, target_start, target_end

            # TODO: Check the sequence with the BLAT result of the query.
            # PSLX format should be used for above check.
            print "QSEQ", query_sequence  # Reverse

            print "TSEQ", target_record.seq
