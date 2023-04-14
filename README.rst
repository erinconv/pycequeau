pycequeau
=======================================

Pycequeau is a library that obtains the necessary files to build and run the CEQUEAU hydrological-water temperature model.

This library was developed and tested in python 3.7, given incompatibilities with GDAL. In order to install it and make it run, we recommend using Anaconda and creating a virtual environment:

.. code-block:: bash
    
    conda create -n pycequeau python=3.7


Now, clone the repository into your local machine and install the dependencies:

.. code-block:: bash
    
    git clone git@github.com:erinconv/pycequeau.git
    cd pycequeau
    conda install --file requirements.txt