# -*- coding: utf-8 -*-
# Generated by Django 1.9.3 on 2016-03-10 12:35
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('grna', '0001_initial'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='guiderna',
            unique_together=set([('target_hit', 'pam_seq', 'cut_position', 'spacer')]),
        ),
    ]