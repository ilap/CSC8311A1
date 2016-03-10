import sys
import tempfile
import regex as re

from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from subprocess import Popen
from subprocess import PIPE

from .models import *


class GuideRNAManager:

    # Global private attributes
    # The default guide RNA spacer length.
    _spacer_length = 20

    _initialised = False

    # Instances of the selected model objects
    _nuclease = None
    _species = None
    _target = None
    _pam = None

    # New models for the objects created by BLAT and search gRNA.
    _target_hit_model = None
    _guide_rna_model = None

    # List of the target gene/sequence found in the genome.
    _hits = None

    # Stream offset for visualising the found gRNA
    _stream_offset = 10

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
            raise Exception('Platform {} is not currently supported.'
                            'Exiting.'.format(self._arch))

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

        # DEBUG: Blat result
        # print blat_qresult

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
            # DEBUG
            # print "LENGTH", float(length), hit[0].hit_start,
            # hit[0].hit_end

            # Strand on HSPFragment, meaning: 0: 5-3, 1: 3-5
            strand = hit[0][0].hit_strand

            # hsp = hit[0]
            # hsp.match_num/length
            hit_score = hit[0].match_num/float(length)

            upstream = int(self._target.upstream)
            downstream = int(self._target.downstream)
            genome_length = int(self._species.length)
            # DEBUG print "SCORE", hit_score

            # TODO: check the genome boundaries.
            # If hit start smaller then upstream then the start must be 0.
            if hit[0].hit_start < upstream:
                target_position = 0
            else:
                target_position = hit[0].hit_start - upstream

            # If hit_end + downstream larger than size of genom the length
            # should be the genome length - hit start
            if (hit[0].hit_end + downstream) > genome_length:
                target_length = genome_length - target_position
            else:
                target_length = length + downstream + downstream

            hit_model = TargetHit.objects.create(species=self._species,
                                                 target=self._target,
                                                 position=target_position,
                                                 strand=strand,
                                                 score=hit_score,
                                                 length=target_length)

            hit_model.save()

            return True

    def initialise_run(self, request):

        if self._initialised:
            raise Exception("The GuideRNABioutisl should be singleton")

        species_id = request.POST['species']
        target = request.POST['target']
        pam_id = int(request.POST['pam'])
        upstream = request.POST['upstream']
        downstream = request.POST['downstream']
        is_nickase = request.POST['is_nickase']
        nuclease_id = request.POST['nuclease']
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
        self._nuclease = Nuclease.objects.get(pk=nuclease_id)

        # print "PAM ID", pam_id, type(pam_id)
        if pam_id < 0:
            self._pam = PAM.objects.filter(nuclease=nuclease_id)
        else:
            self._pam = PAM.objects.filter(nuclease=nuclease_id, pk=pam_id)

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
        # DEBUG print "QUERY", query_file
        # DEBUG print "TARGET", target_file

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
            # print "HIT STRAND", hit.strand
            if hit.strand == 1:
                first_strand = 'plus'
                second_strand = 'minus'
            else:
                first_strand = 'minus'
                second_strand = 'plus'

            hit_start = hit.position
            hit_end = hit.position+hit.length

            # print "HINT", hit, hit_start, hit_end, hit.length,
            # hit_end-hit_start

            # TODO: Check the sequence with the BLAT result of the query.
            # PSLX format should be used for above check.

            # The sequenc is extended by upstream and downstream.
            # query_sequence = query_record.seq[hit_start:hit_end]

            pam_patterns = self._build_pams(self._pam, hit.strand)

            # Find PAMs in + strand
            self._find_pams(hit, query_record.seq, hit_start, hit_end,
                            first_strand, pam_patterns)

            # Find PAMs in - strand
            self._find_pams(hit, query_record.seq, hit_start, hit_end,
                            second_strand, pam_patterns)

            # Result the TargetHits and the GuideRNAs
            target_hits = TargetHit.objects.filter(species=self._species,
                                                            target=self._target)
            return target_hits

    def _build_pams(self, pams, strand):
        # The result is an dict with pams in regex for each strand e.g.
        result = {}

        plus_array = [str(pam).replace('N', '.') for pam in pams]
        minus_array = [str(Seq(str(pam), generic_dna).reverse_complement()
                           ).replace('N', '.') for pam in pams]

        plus_pattern = '(?=(' + '|'.join(plus_array) + '))'
        minus_pattern = '(?=(' + '|'.join(minus_array) + '))'

        if strand == 1:
            result['plus'] = plus_pattern
            result['minus'] = minus_pattern
        else:
            # TODO: Check and unit test this.
            result['minus'] = plus_pattern
            result['plus'] = minus_pattern

        # DEBUG print "DICT", result

        return result

    # TODO: Currently, only supports the same length PAMs. e.g NAG ANG NGG etc.
    def _find_pams(self, hit, query_sequence, start, end, strand, patterns):

        PAM_LENGTH = 3  # TODO: Should be the length of PAM

        if strand == 'plus':
            spacer_start = -self._spacer_length
            spacer_end = 0

        else:
            spacer_start = PAM_LENGTH
            spacer_end = PAM_LENGTH + self._spacer_length

        print "PATTERNS", str(patterns[strand])
        pam_pattern = re.compile(patterns[strand])

        for match in re.finditer(pam_pattern, str(query_sequence[start:end]),
                                 overlapped=True):

            pos = match.start(0)
            pam = match.group(1)

            grna_start = start + pos+spacer_start
            grna_end = start + pos+spacer_end

            # 1. PAM
            ###
            # For comparison
            found_pam_seq = match.group(1)

            # TODO: make PAM length variable.
            temp_pam_seq = match.group(1)
            pam_seq = query_sequence[start + pos:start + pos + PAM_LENGTH]
            pam_seq_c = str(Seq(pam, generic_dna).complement())

            if temp_pam_seq != pam:
                assert Exception("PAM not the same.")

            # 2. SPacer/Target sequence
            ###
            test_spacer_seq = query_sequence[grna_start:grna_end]
            spacer = str(test_spacer_seq)
            spacer_c = str(Seq(spacer, generic_dna).complement())

            # 3. Up/Down stream
            ###
            if strand == 'plus':
                # cut position is always calculated from the sense strand
                cut_position = start + pos - self._nuclease.cut_offset
                up_seq_end = grna_start
                down_seq_start = start + pos+PAM_LENGTH
            else:
                cut_position = start + pos+self._nuclease.cut_offset
                up_seq_end = start + pos
                down_seq_start = grna_end

            up_seq_start = up_seq_end-self._stream_offset
            down_seq_end = down_seq_start + self._stream_offset

            # _c means complement
            up_seq = str(query_sequence[up_seq_start:up_seq_end])
            up_seq_c = str(Seq(up_seq, generic_dna).complement())

            down_seq = str(query_sequence[down_seq_start:down_seq_end])
            down_seq_c = str(Seq(down_seq, generic_dna).complement())

            # DEBUG
            # if strand == 'plus':
            #    print "POS,SEQ", pam_seq, query_sequence[
            #                          start+pos-30:start+pos+3+10]
            # else:
            #    print "POS,SEQ", pam_seq, query_sequence[
            #                          start+pos-10:start+pos+3+20+10]

            grna = GuideRNA.objects.create(nuclease=self._nuclease,
                                           target_hit=hit,
                                           pam_seq=pam_seq,
                                           pam_seq_c=pam_seq_c,
                                           spacer=spacer,
                                           spacer_c=spacer_c,
                                           up_seq=up_seq,
                                           up_seq_c=up_seq_c,
                                           down_seq=down_seq,
                                           down_seq_c=down_seq_c,
                                           cut_position=cut_position,
                                           is_sense_strand=int(strand ==
                                                               'plus'))

            grna.save()

        '''
        # DEBUG print "GUIDERNA", "STRAND", pam_seq, strand
        if strand == 'plus':
            print grna.target_hit
            # Debug found sequence
            print "+ Seq:"
            print grna.up_seq, grna.spacer, grna.pam_seq, grna.down_seq

            print "- Seq:"
            print grna.up_seq_c, grna.spacer_c, grna.pam_seq_c, grna.down_seq_r
            print "POS", start+pos, grna.cut_position
            print "Is sense/+", grna.is_sense_strand

            print "SEQUENCE"
            print query_sequence[grna.cut_position-26:grna.cut_position+17]
        else:
            print grna.target_hit
            # Debug found sequence
            print "+ Seq:"
            print grna.up_seq, grna.pam_seq, grna.spacer,  grna.down_seq

            print "- Seq:"
            print grna.up_seq_c, grna.pam_seq_c, grna.spacer_c,  grna.down_seq_r
            print "POS", start+pos, grna.cut_position
            print "Is sense/+", grna.is_sense_strand

            print "SEQUENCE"
            print query_sequence[grna.cut_position-14:grna.cut_position+29]
            '''