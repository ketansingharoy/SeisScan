============
Installation
============

SeisScan is currently running on Mac OS. SeisScan runs on Python 3.10 and up. We recommend you use the latest version of python 3 if possible.



Dependencies
============
* python [>= 3.10]
* obspy [>=1.4.0]
* pandas [>=2.2.2]
* numpy [>=1.26.4]
* utm [>=0.7.0]
* dask [>=2024.8.0]
* distributed [>=2024.8.0]
* jupyter [>=1.0.0]

Install via Pip
===============
.. code-block:: bash

    $ pip install seisscan


Install via Anaconda
====================
.. code-block:: bash

    $ conda create -n seisscan-env python=3.10
    $ conda activate seisscan-env
    $ conda install -c ksr22 seisscan