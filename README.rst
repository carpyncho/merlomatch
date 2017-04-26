MMatch
======

A Cross-Match functionality for VVV sources


Features
--------

*   Command line interface with standard utilities
*   Multiprocess



Support
-------

-   **Contact:** jbc.develop@gmail.com
    https://groups.google.com/forum/#!forum/corral-users-forum
-   **Issues:** If you have issues please report them as a issue
    here: https://github.com/carpyncho/mmatch/issues


User Installation
-----------------

The easiest way to install mmatch is using pip

.. code-block:: bash

    $ pip install -U mmatch


Developmment
------------

Install project by cloning from `mmatch github <https://github.com/carpyncho/mmatch/>`__:

.. code-block:: bash

    $ git clone https://github.com/carpyncho/mmatch.git

and by making ``pip install -e .``, or the classic ``python setup.py install``.

You can also run the install by giving the github link directly to pip:

.. code-block:: bash

    $ pip install -e git+https://github.com/carpyncho/mmatch.git


Help
----

.. code-block:: bash

    usage: mmatch.py [-h] [-v]
                 (--sources-file PATH | --sources "RA DEC" ["RA DEC" ...])
                 [--quiet] --pawprints PATH [PATH ...] [--radius RADIUS]
                 --vvv-flx2mag PATH [--working-directory PATH] [--output PATH]
                 [--procs NUMBER]

    A Cross-Match functionality for VVV sources

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      --sources-file PATH, -srcf PATH
                            CSV file with at least 2 columns "ra, dec"
      --sources "RA DEC" ["RA DEC" ...], -srcs "RA DEC" ["RA DEC" ...]
                            Quoted ra and dec parameters example: "17:29:21.4
                            -30:56:02"
      --quiet, -q           set log level to warning
      --pawprints PATH [PATH ...], -pwps PATH [PATH ...]
                            Pawprint files or directory
      --radius RADIUS, -r RADIUS
                            Radiuous to make the crossmatch
      --vvv-flx2mag PATH, -f2m PATH
                            vvv_flx2mag command PATH
      --working-directory PATH, -wd PATH
                            directory to store the temporary files
      --output PATH, -o PATH
                            destination of your pawprint
      --procs NUMBER, -p NUMBER
                            Number of processors to use (0 run in the same
                            process)

    BSD-3 Licensed - IATE-OAC: http://iate.oac.uncor.edu/



License
-------

-   BSD-3: https://opensource.org/licenses/BSD-3-Clause
-   BSD-3 Licence Explained: https://tldrlegal.com/license/bsd-3-clause-license-(revised)
