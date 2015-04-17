# Copyright (C) 2014, 2015  Peter Batzel and Jason Sydes
#
# This file is part of Prost!.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the License with this program.
#
# Written by Jason Sydes and Peter Batzel.
# Advised by Thomas Desvignes.

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

# peewee - python ORM
from peewee import *
from playhouse.proxy import Proxy

# For peewee datetimes
import datetime

# The peewee database proxy
DB_PROXY = Proxy()


###############
### Classes ###
###############

class SequenceDB(Model):
    """Database table to cache the sequence database used in the genomic
    alignment.

    For example, if using BBMap, this will be a BBMap database.
    """

    name = CharField()
    created_at = DateTimeField(default = datetime.datetime.now, null = False)

    class Meta:
        database = DB_PROXY
        db_table = 'sequence_dbs'


class ShortSeqCache(Model):
    """Database table to cache aspects of ShortSeq objects.

    This table is designed to store the results of the genomic alignment as
    well as the results of the designation step one.
    """
    db = ForeignKeyField(SequenceDB,
            db_column = 'sequence_db_id',
            related_name='short_seqs')
    sequence = CharField(
            unique = True,
            index = True,
            null = False)
            # max_length = COL_LEN__SHORT_SEQS__SEQUENCE)
    designation_integer = IntegerField(null = False)
    created_at = DateTimeField(default = datetime.datetime.now, null = False)

    class Meta:
        database = DB_PROXY
        db_table = 'short_seqs'


class GenomicLocationCache(Model):
    """
    Database table to cache the genomic locations of ShortSeq objects.
    """
    short_seq = ForeignKeyField(ShortSeqCache,
            db_column = 'short_seq_id',
            related_name='genomic_locations')
    lg = CharField(null = False)
    start = IntegerField(null = False)
    end = IntegerField(null = False)
    # The designation stores true designtations for 3 and above only.
    designation = CharField(
            null = True)
            # max_length = COL_LEN__GENOMIC_LOCATION__DESIGNATION)
    created_at = DateTimeField(default = datetime.datetime.now, null = False)

    class Meta:
        database = DB_PROXY
        db_table = 'genomic_locations'
        indexes = (
            # create a unique index on short_seq_id/lg/start/end
            (('short_seq', 'lg', 'start', 'end'), True),
        )


# vim: softtabstop=4:shiftwidth=4:expandtab
