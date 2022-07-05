# genIAL
is a python code that generates intrinsic alignment (IA) of galaxies by assigning intrinsic shapes and orientations based on the galaxies' photometric properties as well as the orientation and angular momentum of their dark matter host haloes. The code is fully vectorized, allowing it to be run on large datasets using the Apache Spark framework.

The modeling is done in three steps:

    - assigning 3D galaxy axis ratios
    - assigning 3D galaxy orientations
    - projecting 3D galaxies along the observed line of sight to obtain the 2D ellipticities

# reference
Details about the IA model are described here https://arxiv.org/abs/2206.14219.
Please cite this paper when using the code for your research.

# online example
An online demo notebook showing how to run genIAL on an input file from a cosmological simulation can be found here:
https://colab.research.google.com/drive/1NUVA05EiyoGfrNyxVAS_2bZSOPVq5ngZ?usp=sharing

# contact
Kai Hoffmann: kai (dot) d (dot) hoffmann (at) gmail (dot) com
