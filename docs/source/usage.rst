Usage
=====

.. _installation:

Installation
------------

To use ReadOIFITS, first install it using conda:

.. code-block:: console

   conda install readoifits

Reading your oifits files
-------------------------

To read your oifits files,
you can use the ``ReadOIFITS.read()`` function:

.. autofunction:: ReadOIFITS.read()

The ``directory``  should be a string with the path to the directory with the
data files and ``files`` should be a string with the name of the files.
It is possible to specify ``*fits`` to select all the files ending with fits
in the specified directory.

For example:

>>> import ReadOIFITS as oifits
>>> data_dir = '/Users/XXX/path/to/you/data/'
>>> files = '*.fits'
>>> oifits.read(data_dir, files)
