from __future__ import unicode_literals
from django.db import models

# Create your models here.

class GuideRNA(models.Model):
	species = models.CharField(max_length=200)

	target = models.CharField(max_length=2000)


	upstream = models.CharField(max_length=10)
	downtream = models.CharField(max_length=10)

	pam = models.CharField(max_length=10)
	cas_type = models.CharField(max_length=10)

	def __str__(self):
		return self.species



