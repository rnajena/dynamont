===============================
Dynamont
===============================

.. image:: https://img.shields.io/travis/JannesSP/dynamont.svg
        :target: https://travis-ci.org/JannesSP/dynamont
.. image:: https://circleci.com/gh/JannesSP/dynamont.svg?style=svg
    :target: https://circleci.com/gh/JannesSP/dynamont
.. image:: https://codecov.io/gh/JannesSP/dynamont/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/JannesSP/dynamont

Segmentation/resquiggling tool for ONT signals.

Installation
------------

Pypi/pip
~~~~~~~~

.. code:: bash

   pip install dynamont

Conda
~~~~~

.. code:: bash

   conda install mamba
   mamba create -n dynamont -c jannessp dynamont
   conda activate dynamont

--------------

Usage
-----

.. code:: bash
   # segment a dataset
   dynamont-resquiggle -r <path/to/pod5/dataset/> -b <basecalls.bam> --mode basic --model_path <path/to/model> -o <output.csv> -p <pore>

   # train model
   dynamont-train -r <path/to/pod5/dataset/> -b <basecalls.bam> --mode basic --model_path <path/to/init/model> -o <output/path> -p <pore>

-----