.. image:: https://readthedocs.org/projects/pysoftwhere/badge/?version=latest

pySoftWhere:  Automated Interface Analysis of Soft Matter Nanostructures
=========================================================================

**pySoftWhere** is a Python library for analyzing interfaces of soft matter nanostructures, such as *micelles, nanoparticles, monolayers, and polymer films*. N.B.: This project is still in development!	

Using pySoftWhere
-----------------

Please consult the documentation on `Read the Docs <https://pysoftwhere.readthedocs.io/en/latest/index.html>`_

Installation
------------

1. It is recommended to create a new virtual environment in which the pySoftWhere dependencies can be installed, for example

.. code-block:: console
   
   conda create --name pysw_env
   conda activate pysw_env


2. Make a directory to download pySoftWhere and change to this directory

.. code-block:: console
   
   (pysw_env) mkdir pysw_dir
   (pysw_env) cd pysw_dir

3. Download pySoftWhere from Github

.. code-block:: console
   
   (pysw_env) git clone https://github.com/rmziolek/pySoftWhere.git 

4. Install pySoftWhere and its dependencies as required

.. code-block:: console

   (pysw_env) cd pySoftWhere
   (pysw_env) pip install -e .

5. You may need to add the pySoftWhere installation directory to your .bashrc or .zshrc file

.. code-block:: console
    
    PATH=$PATH:/<path to>/pysw_dir	

That's it, you're ready to use pySoftWhere (remember to activate pysw_env!)
