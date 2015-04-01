=========
pyrefflat
=========

Pyrefflat is a parser for refFlat files, used by UCSC to annotate genes and their exons, in pure python.
It includes a reader and writer object for easy manipulation of refFlat files.

Installation
------------

To install pyrefflat, clone the repository at  https://git.lumc.nl/a.h.b.bollen/pyrefflat .
Run ``python setup.py install`` to install the module. It is recommended you use a virtual environment.

Usage
-----

Reading refFlat files
~~~~~~~~~~~~~~~~~~~~~
Pyrefflat provides a ``Reader`` object for reading refFLat files.

.. code-block:: python

    from pyrefflat import Reader
    reader = Reader(filename="file.reFlat")

A reader is an iterator that returns records. Each record has associated exons.
E.g., to print the start site of every exon in every record, one would do:

.. code-block:: python

    for record in reader:
        for exon in record.exons:
            print exon.start


Writing refFlat files
~~~~~~~~~~~~~~~~~~~~~
Writing refFlat files is done with the ``Writer`` object, which consumes ``Record`` instances.
E.g., to copy a refFlat file using a ``Reader`` and ``Writer`` one could do:

.. code-block:: python

    from pyrefflat import Reader, Writer
    reader = Reader(filename="original.refFlat")
    writer = Writer(filename="writer.refFlat")

    for record in reader:
        writer.write(record)
    writer.close()
    reader.close()


**NOTE: this part is still in flux, so liable to change!**

New records can be generated using the ``RecordFactory`` object in ``pyrefflat.factories``.
Initializing a ``RecordFactory`` object without any arguments will result in an empty line.
If one initializes a ``RecordFactory`` with a ``Record`` object, it will use the ``Record`` as a template.
Fields can be set by using the a setter method. Fields can be imported from ``pyrefflat.generics``.
The following creates a record on chromosome 1, with gene ``AAAA``, ranging from position 1 to 1000 and having 1 exon.

.. code-block:: python

    from pyrefflat.factories import RecordFactory
    from pyrefflat.generics import NUMERIC_COLUMNS, NUMERIC_LIST_COLUMNS, STRING_COLUMNS

    factory = RecordFactory()
    factory.set_gene("AAAA")
    factory.set_transcript("AAAA")
    factory.set_chromosome(1)
    factory.set_strand( "+")
    factory.set_transcription_start(1)
    factory.set_transcription_end(1000)
    factory.set_cds_start(1)
    factory.set_cds_end(1000)
    factory.set_exon_count(1)
    factory.set_exon_starts([1])
    factory.set_exon_ends([1000])

    record = factory.make()


Tools
-----

Apart from a parser, it includes several standalone tools. These additionally depend on the ``pysam`` and ``vcf`` modules.

gVCFCoverage
~~~~~~~~~~~~
This tool calculates coverage on the regions specified in the refFlat exon fields in a gVCF file.
It can output in three different formats, namely simple tab-delimited (csv or tsv) format, BED format or as a JSON.
It supports using the GQ field, allowing to filter only those regions with a minimum GQX value.

createMargin
~~~~~~~~~~~~
This tool adds a margin around each exon and writes the result to a new refFlat file.

refFlat2Bed
~~~~~~~~~~~
This tool converts a refFlat file to a BED file, with the regions based on the exons.


License
-------
pyrefflat is MIT licensed.