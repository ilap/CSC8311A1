from django.test import TestCase

from .models import Nuclease
from .models import PAM

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


class PAMTestCase(TestCase):

    def setUp(self):
        self._nuclease = Nuclease.objects.create(name="Cas9",
                                                 is_nickase=False,
                                                 is_downstream=False,
                                                 cut_offset=4)

        PAM.objects.create(nuclease=self._nuclease, pam="NGG")
        PAM.objects.create(nuclease=self._nuclease, pam="NAG")

        self._sequence = Seq("AGTACAGAGCTGGTACGGTTTGGGC", generic_dna)

    def tearDown(self):
        pass

    def test_sequence(self):
        # Test regex
        self.assertEqual(self._sequence.find("GG"), 11)

    def test_pam_reversecomplements(self):

        pam = PAM.objects.get(pam="NGG")
        self.assertEqual(Seq(str(pam.pam), generic_dna).reverse_complement(),
                         "CCN")

        pam = PAM.objects.get(pam="NAG")
        self.assertEqual(Seq(str(pam.pam), generic_dna).reverse_complement(),
                         "CTN")

    def test_pam_simple_regex(self):

        pams = PAM.objects.all()

        patterns = [str(pam.pam) for pam in pams]
        plus_pattern = "(?=(" + '|'.join(patterns).replace('N', '.') + "))"

        patterns = []
        for pam in pams:
            tmp_str = str(Seq(pam.pam).reverse_complement()).replace("N", ".")
            patterns.append(tmp_str)

        minus_pattern = "(?=(" + '|'.join(patterns).replace('N', '.') + "))"

        self.assertEqual(plus_pattern, "(?=(.GG|.AG))")
        self.assertEqual(minus_pattern, "(?=(CC.|CT.))")

        # print "+ PATTERN", plus_pattern
        # print "- PATTERN", minus_pattern

        # pattern = re.compile(plus_pattern)

        # for alma in re.finditer(pattern, str(self._sequence),
        # overlapped=True):
        # print "AAAA", alma.group(1), alma.start(0), alma.end(0)
