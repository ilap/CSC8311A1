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

    # Target position in species if exists
    position = models.IntegerField(validators=[MinValueValidator(-1)])

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
        self.position = -1

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

    def get_fasta_file(self):
        return self.sequence_file.path


# Genome of a species to find gRNAs by target sequence(s)
class Species(models.Model):

    name = models.CharField(max_length=100, unique=True)
    fasta_file = models.FileField(storage=gfs)
    genbank_file = models.FileField(storage=gfs)
    # Sequence length
    length = models.PositiveIntegerField()

    def __str__(self):
        return self.name

    def get_fasta_file(self):
        return self.fasta_file.path

    def get_genbank_file(self):
        return self.genbank_file.path


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


# Model for gRNAs
class GuideRNA(models.Model):
    species = models.ForeignKey(Species, null=True)
    target = models.ForeignKey(Target, null=True)
    pam = models.ForeignKey(PAM, null=True)
    nuclease = models.ForeignKey(Nuclease, null=True)

    grna = models.CharField(max_length=20, null=True)

    # Position of found gRNA in Genome and/or Target sequence.
    genome_position = models.BigIntegerField(default=-1, null=True)
    target_position = models.BigIntegerField(default=-1, null=True)

    class Meta:
        unique_together = (("species", "nuclease", "target",
                            "pam"),)

    def __str__(self):
        return self.grna
