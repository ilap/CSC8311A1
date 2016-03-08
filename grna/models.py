from __future__ import unicode_literals
from django.db import models
from django.core.files.storage import FileSystemStorage
from django.core.validators import MinValueValidator
from django.core.files.base import ContentFile
from Bio import SeqIO

import hashlib

tfs = FileSystemStorage(location='grna/targets')
gfs = FileSystemStorage(location='grna/genomes')


# Create your models here.
# Target sequence to find gRNA in a genome of a species
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


# Genome of a species to find gRNAs by target sequence(s)
class Species(models.Model):

    name = models.CharField(max_length=100, unique=True)
    fasta_file = models.FileField(storage=gfs)
    genbank_file = models.FileField(storage=gfs)
    # Sequence length
    length = models.PositiveIntegerField()

    def __str__(self):
        return self.name


# Model for CRISPR/Cas9 or CRISPR/Cpf1 endonucleases
class Nuclease(models.Model):
    name = models.CharField(max_length=10, default='Cas9')
    is_nickase = models.BooleanField(default=True)

    is_downstream = models.BooleanField(default=True)
    cut_offset = models.PositiveIntegerField(validators=[MinValueValidator(
        1)], null=True)

    def __str__(self):
        return self.name


# Model for PAM (Protospacer Adjacent Motif)
class PAM(models.Model):
    nuclease = models.ForeignKey(Nuclease)
    pam = models.CharField(max_length=20)

    class Meta:
        unique_together = (("nuclease", "pam"),)

    def __str__(self):
        return self.pam


# Model for Target Hits in the genome
class TargetHit(models.Model):
    species = models.ForeignKey(Species, null=True)
    target = models.ForeignKey(Target, null=True)

    length = models.BigIntegerField(default=-1, null=True)

    # Target position in genome of species if exists
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


# Model for gRNAs
class GuideRNA(models.Model):
    target_hit = models.ForeignKey(TargetHit, null=True)
    pam = models.ForeignKey(PAM, null=True)

    grna = models.CharField(max_length=20, null=True)

    # Position of found gRNA in Genome
    # gRNA position is Target is calculated from TargetHit's position
    # and gRNA position in Genmee.
    position = models.BigIntegerField(default=-1, null=True)

    class Meta:
        unique_together = (("target_hit", "pam", "position"),)

    def __str__(self):
        return self.grna
