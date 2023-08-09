Installing pySoftWhere
=====


pySoftWhere is available to download from `GitHub <https://github.com/rmziolek/pySoftWhere>`_

1. It is recommended to create a new virtual environment in which the pySoftWhere dependencies can be installed, for example

.. code-block:: console
   
   conda create --name pysw
   conda activate pysw


2. Make a directory to download pySoftWhere and change to this directory

.. code-block:: console
   
   (pysw) mkdir pysoftwhere
   (pysw) cd pysoftwhere

3. Download pySoftWhere from Github

.. code-block:: console
   
   (pysw) git clone https://github.com/rmziolek/pySoftWhere.git 

4. Install pySoftWhere and its dependencies as required

.. code-block:: console
   
   (pysw) pip install -e .

5. You may need to add the installation directory path to your .bashrc or .zshrc file

.. code-block:: console
    
    PATH=$PATH:/Users/<path to>/pySoftWhere	

That's it, you're ready to use pySoftWhere

