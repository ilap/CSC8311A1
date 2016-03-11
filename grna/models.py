from __future__ import unicode_literals
from django.db import models
from django.core.files.storage import FileSystemStorage
from django.core.validators import MinValueValidator
from django.core.files.base import ContentFile
from Bio import SeqIO

import hashlib

tfs = FileSystemStorage(location='grna/targets')
gfs = FileSystemStorage(location='grna/genomes')



# Target is a model for target sequences or gene sequences for searching guide
# RNAs in similar or homologous sequence in a genome.
# One genome can have more than one genes or homologoues sequences and the
# TargetHit models (see further below) is the model for those founds/hits.
class Target(models.Model):

    # TODO: handle the unique targets.
    name = models.CharField(max_length=100)

    # Target sequence in FASTA, e.g. gene, manually uploaded sequence, or
    # sequence from offsets of genome
    sequence_file = models.FileField(storage=tfs)

    # Target length
    length = models.BigIntegerField(default=-1, null=True)

    # If the Target is found in the genome of species, then it specifies
    # the extension of the sequence to find gRNA in the extended target
    # by the up/down stream.
    upstream = models.PositiveIntegerField(validators=[MinValueValidator(0)])
    downstream = models.PositiveIntegerField(validators=[MinValueValidator(0)])

    def __str__(self):
        return self.name

    def init_target(self, target, upstream, downstream):
        self.upstream = upstream
        self.downstream = downstream

        # The file name would be the hash of the sequence. MD5 is suitable as
        # it is not for encrypting but for unique name for a file based on
        # the sequene
        h = hashlib.md5()
        h.update(target)
        hash_name = h.hexdigest() + ".fa"

        # TODO: Handle invalid contents, means no fasta format.
        self.sequence_file.save(hash_name, ContentFile(target))

        # TODO: Handle no fasta format file and more than one fasta records
        # currently it handles only one FASTA record

        # Get the id from fasta file and set the model name into it
        fasta_file = self.sequence_file.path
        seq_record = SeqIO.read(fasta_file, "fasta")

        if seq_record:
            self.name = seq_record.id
        else:
            self.name = hash_name

        self.save()


# Species is the model of the species genom where the guide RNAs will be
# searched in the Hits (see detail in TargetHit model further below) in the
# chose genome of  species.
class Species(models.Model):

    name = models.CharField(max_length=100, unique=True)
    fasta_file = models.FileField(storage=gfs)
    genbank_file = models.FileField(storage=gfs)
    # Sequence length
    length = models.PositiveIntegerField()

    def __str__(self):
        return self.name


# Model for CRISPR/Cas9 or CRISPR/Cpf1 endonucleases which binds to the
# specific sequence (PAM) in the genome and cleavage the DNA.
class Nuclease(models.Model):
    name = models.CharField(max_length=10, default='Cas9')

    # Nickase means create nicks instead of a blunt cut of the fully
    # functional Cas9 enzyme.
    is_nickase = models.BooleanField(default=True)

    # The spacer sequence of the guide RNA relative to the PAM sequence it
    # can be either right up or down next to PAM sequence.
    is_downstream = models.BooleanField(default=True)

    # Where the enzyme cuts the DNA, it can be upstream or downstream
    # relative to the PAM sequence
    cut_offset = models.PositiveIntegerField(validators=[MinValueValidator(
        -10)], null=True)

    def __str__(self):
        return self.name


# Model for PAM (Protospacer Adjacent Motif) sequence where the Cas9/Cpf1
# endonuclease binds and cleavage the DNA
class PAM(models.Model):
    nuclease = models.ForeignKey(Nuclease)

    # The regex of the PAM sequence of the specified nuclease for example
    # NGG, NAG, or NNGNNGG, whihc depends on the enzyme.
    pam = models.CharField(max_length=20)

    class Meta:
        unique_together = (("nuclease", "pam"),)

    def __str__(self):
        return self.pam


# Model for Target Hits
# Target hits the sequences similar to the target (gene or any selected)
# sequence in the query sequence (genome).
# More than one homologous of the target can be found in the queried genome,
# therefore the guide RNAs must be found in these finds (hits)
# In default the homology is set to 90% in the BLAT whihc finds homologous
# sequences (here the target sequences) in a query sequence (here the species
# genome)
class TargetHit(models.Model):
    # The species object where the target is found.
    species = models.ForeignKey(Species, null=True)

    # The original target object
    target = models.ForeignKey(Target, null=True)

    length = models.BigIntegerField(default=-1, null=True)

    # Target position in genome if exists
    position = models.BigIntegerField(default=-1, null=True, validators=[
        MinValueValidator(-1)])

    # Strand can be either plus, sense or coding strand (means 5' to 3' or
    # minus, antisense or not coding strand
    #
    ORIENTATION = (
        (1, 'Coding Strand'),
        (-1, 'Template Strand')
    )
    strand = models.IntegerField(default=1, choices=ORIENTATION, null=True)

    # Hit Score would be larger than 90% as BLAT default value for search
    # TODO: Parameterise this property in start page
    score = models.FloatField(default=0.9, null=True)

    class Meta:
        unique_together = (("species", "target", "position"),)


# Model for Guid RNAs found in similar sequence to the target gene/sequence
#  in the genome of specified species.
class GuideRNA(models.Model):
    # The target hit (the homologous gene or other sequence(s) of the queried
    # genome/sequence).
    target_hit = models.ForeignKey(TargetHit, null=True)

    # The used endonuclease for generating guide RNAs
    nuclease = models.ForeignKey(Nuclease, null=True)

    # PAM sequence of the guide RNA
    pam_seq  = models.CharField(max_length=20, null=True)
    pam_seq_c = models.CharField(max_length=20, null=True)

    # Guide RNA Target/Spacer Sequence, the nuclease cut the DNA based on
    # this sequences.
    spacer = models.CharField(max_length=20, null=True)
    spacer_c = models.CharField(max_length=20, null=True)

    # Up/down stream sequence relative to the guide RNA (PAM + spacer)
    up_seq = models.CharField(max_length=20, null=True)
    up_seq_c = models.CharField(max_length=20, null=True)
    down_seq = models.CharField(max_length=20, null=True)
    down_seq_c = models.CharField(max_length=20, null=True)

    # Position of found gRNA in Genome
    # gRNA position is Target is calculated from TargetHit's position
    # and gRNA position in Genome.
    cut_position = models.BigIntegerField(default=-1, null=True)

    # is the Guide RNA on sense/plus strand?
    is_sense_strand = models.BooleanField(default=True)

    class Meta:
        unique_together = (("target_hit", "pam_seq", "cut_position", "spacer"),)

    def __str__(self):
        return self.spacer

