# full_pysoftk


This is the full version of pysoftk that includes:

1. ring ring and all other functions
2. umap and hbdscan
3. intrinsic densities (linear and spherical)


# How to install it

you need to install first pysoftwhere and then afterwards pysoftk. Using either conda or virtual enviorements, you can use **pip install -e .**
in both cases. In the case of conda, please enable the environment. For a virtual enviorement, please just source the **activate** file.


# To run the test and the pytest inside python

Due to the size of the trajectory file (115MB), I have to split the trajetory into 2 file which are located in the folder **test**. 
To join the trajectory just type (**cat intrinsic_density_movie_split.aa intrinsic_density_movie_split.ab > intrinsic_density_movie.xtc**)
this enables to use the file **fail.py**. 

If you want to run all tests within the folder **../pysoftk/test/test_pol** , please move the trajectory **.xtc** to the folder
**../pysoftk/test/test_pol/data** and then type inside the folder **../pysoftk/test/test_pol** the command **pytest**.
