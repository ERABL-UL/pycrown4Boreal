[![manaakiwhenua-standards](https://github.com/manaakiwhenua/pycrown/workflows/manaakiwhenua-standards/badge.svg)](https://github.com/manaakiwhenua/manaakiwhenua-standards)


# PyCrown - Fast raster-based individual tree segmentation for LiDAR data
Author: Dr Jan Schindler (formerly Zörner) (<mailto:schindlerj@landcareresearch.co.nz>)

Published under GNU GPLv3


# Summary
PyCrown is a Python package for identifying tree top positions in a canopy height model (CHM) and delineating individual tree crowns.

The tree top mapping and crown delineation method (optimized with Cython and Numba), uses local maxima in the canopy height model (CHM) as initial tree locations and identifies the correct tree top positions even in steep terrain by combining a raster-based tree crown delineation approach with information from the digital surface model (DSM) and terrain model (DTM).

*Citation:*

Zörner, J.; Dymond, J.; Shepherd J.; Jolly, B. PyCrown - Fast raster-based individual tree segmentation for LiDAR data. Landcare Research NZ Ltd. 2018, https://doi.org/10.7931/M0SR-DN55

*Research Article:*

Zörner, J., Dymond, J.R., Shepherd, J.D., Wiser, S.K., Bunting, P., Jolly, B. (2018) Lidar-based regional inventory of tall trees - Wellington, New Zealand. Forests 9, 702-71. https://doi.org/10.3390/f9110702

# Changes I have made to make it work on our data, which are 20x20 plots of boreal forest

= = = = = = = = = = == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

= = Adaptation de l'algorithme "PyCrown - Fast raster-based individual tree segmentation for LiDAR data" = =

= = = = = = = = = = = = = = = = == = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

LÉGENDE:

 [+] indique que l'élément a été ajouté
 
 [-] indique que l'élément a été modifié
 
 [=] indique que l'élément était déjà là et N'A PAS été modifié mais vaut la peine d'être expliqué à des fins de compréhension
 
 [def] Le terme qui suit est le nom de la fonction qui exécute les tâches listées

= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

[+] Tout mettre l'algorithme dans une for loop pour avoir la segmentation de chaque parcelle de 20x20 mètres;

[+] Création DTM et DSM pour l'obtention des CHM:

	[+]­­­­­­ L'integrer dans la for loop pour avoir le CHM de chaque parcelle;
 
	[+] Utiliser tous les points pour le DSM;
 
	[+] Utiliser seulement les points du sol pour le DTM;
 
	[+] CHM = DSM - DTM.

[-] Trouver les sommets des arbres dans la parcelle:

	[-] Trouver les maximums pour chaque ws (window screen, qui est fait un noyau qui glisse sur le CHM) de 3x3 pixels dans le CHM (diminué de 5x5 à 3x3 parce que la parcelles n'est que de 20x20 et ce n'est pas tous les sommets qui se faisaient identifier);
 
	[-] La hauteur minimale d'un sommet est de 2 mètres. Auparavant à 20 mètres, ce qui est beaucoup trop. 2 mètres est bon vu que ca trouve les petits arbres aussi;
 
	[=] On trouve le centre de masse de deux sommets si le même sommet est trouvé côte-à-côte plusieurs fois.
	
[-] Utilisation de dalponteCIRC_numba pour la création des couronnes d'arbres:

	[-] Definition: "Crown delineation based on Dalponte and Coomes (2016) and lidR R-package (https://github.com/Jean-Romain/lidR/). In contrast to the moving window growing scheme from the original algorithm this code implements a circular region growing around the tree top which leads to smoother crown patterns and speeds up the calculation by one order of magnitude";
	
 [-] Prend en entrée les paramètres suivant:
 
		[=] CHM: on sé c koi;
  
		[=] Trees: Coord des sommets des arbres trouvés, nx2;
  
		[=] th_tree: Mis à la valeur de 2. Threshold below which a pixel cannot be a tree. Default 2;
  
		[-] th_seed: Valeur était à 0.45, modifié à 0.4
		  "Growing threshold 1. A pixel is added to a region if its height is greater than the tree height multiplied by this value. It should be between 0 and 1. Default 0.45";
    
		[-] th_crown: Valeur était à 0.55, modifié à 0.5
		  "Growing threshold 2. A pixel is added to a region if its height is greater than the current mean height of the region multiplied by this value. It should be between 0 and 1. Default 0.55.";
    
		[=] max_crown: Laissé à la valeur de 10. "Maximum value of the crown diameter of a detected tree (in pixels). Default 10".
		
	[-] En somme, j'ai diminué les thresholds th_seed et th_crown parce que les couronnes des arbres étaient trop petites avec les valeurs par défaut. PyCrown était utilisé pour des arbres feuillus, grands et assez étroits.
		
[-] Enlever les arbres qui sont de moins de 5 mètres de hauteur (paramètre auparavant établi à 20m). Les couronnes sont mises à jour par la suite. (def screen_small_trees)

[=] Changer les couronnes raster en polygones rasters individuels (def crowns_to_polys_raster):

	[=] Utilise l'outil "shapes" dans le package "rasterio.feature" (installé avec pip, ça n'a pas rapport avec PyCrown).
	
[-] "Smooth crown polygons using Dalponte & Coomes (2016) approach: Builds a convex hull around first return points (which lie within the rasterized crowns). Optionally, the trees in the LiDAR point cloud are classified based on the generated convex hull" (def crowns_to_polys_smooth):

	[-] J'ai modifié pour que le fichier de nuage de point se cré toujours. C'EST DANS CETTE FONCTION QUE LES FICHIERS DE NUAGES SONT CRÉÉS;
 
	[=] Converti le nuage de points en points Shapely;
	[-] Converti les polygones rasters en polygones Shapely. J'ai modifié ça un peu parce que les y des polygones étaient multipliés par -1... Le nuage de points correspondait jamais à aucun polygone. Pas vraiment la peine de le mentionné;
 
	[=] Chaque polygone est ajouté au lot de polygones et on regroupe les points du nuage qui appartiennent au même polygone;
 
	[=] Pour chaque polygone, on prend les points qu'il contient pour créer un nouveau polygone qui suit le contour de ces points. Le polygone est donc lissé et ça cré un masque pour l'attribution d'un index à chaque arbre segmenté;
 
	[-] Pour l'écriture comme tel des fichiers, je l'ai laissé en .las parce que Reza les voulait en .ply. J'ai ajouté les points de sols à chaque nuage de points final segmenté pour pas juste avoir les arbres.




# Purpose and methods
A number of open-source tools to identify tree top locations and delineate tree crowns already exist. The purpose of this package is to provide a fast and flexible Python-based implementation which builds on top of already well-established algorithms.

Tree tops are identified in the first iteration through local maxima in the smoothed CHM.

We re-implement the crown delineation algorithms from **Dalponte and Coomes (2016)** in Python. The original code was published as R-package *itcSegment* (<https://cran.r-project.org/package=itcSegment>) and was further optimized for speed in the *lidR* R-package (<https://cran.r-project.org/package=lidR>).

Our Cython and Numba implementations of the original algorithm provide a significant speed-up compared to *itcSegment* and a moderate improvement over the version available in the *lidR* package.

We also adapted the crown algorithm slightly to grow in circular fashion around the tree top which gives crown a smoother, more natural looking shape.

We add an additional step to correct for erroneous tree top locations on steep slopes by taking either the high point from the surface model or the centre of mass of the tree crown as new tree top.

Reference:

**Dalponte, M. and Coomes, D.A. (2016)** *Tree-centric mapping of forest carbon density from airborne laser scanning and hyperspectral data*. Methods in Ecology and Evolution, 7, 1236-1245.


# Main outputs
* **Tree top locations** (stored as 3D ESRI .shp-file)
* **Tree crowns** (stored as 2D ESRI .shp-file)
* **Individual tree classification of the 3D point cloud** (stored as .las-file)


# Contributors
* Dr Jan Zörner (Manaaki Whenua - Landcare Research, Lincoln, New Zealand)
* Dr John Dymond (Manaaki Whenua - Landcare Research, Palmerston North, New Zealand)
* Dr James Shepherd (Manaaki Whenua - Landcare Research, Palmerston North, New Zealand)
* Dr Ben Jolly (Manaaki Whenua - Landcare Research, Palmerston North, New Zealand)


# Requirements
It is assumed that you generated a canopy height model (CHM), digital surface model (DSM) and digital terrain model (DTM) from the LiDAR dataset before running *PyCrown*.
If you want to classify individual trees in the point cloud, it is recommended to normalize heights to *height above ground elevation* (also done externally).

For processing laser scanning data we recommend the open-source software *SPDLib* (http://www.spdlib.org).


# Installation and environment set-up
**Python 3.6 is required.**

Tested on: Windows 10, Debian 9 (Stretch), Fedora 28, Ubuntu 18.04 & 16.04

## Environment set-up
### With Conda package manager (recommended)
#### Create the environment and install all required packages

`conda env create`

#### Activate the environment

Windows: `activate pycrown-env`

Linux: `source activate pycrown-env`

### With Python's venv and pip
#### Create the environment

`python -m venv pycrown-env`

Linux: `source pycrown-env/bin/activate`

Windows: `pycrown-env\Scripts\activate.bat`

#### Install all required packages

`python -m pip install --upgrade pip`

`pip install -r requirements.txt`

## Run Tests
There are only some rudimentary tests provided at the moment, but it is advised to check that everything works:

`python setup.py test`

## Install PyCrown
Build and install the PyCrown module with:

`python setup.py install`


# Common problems
## laspy.util.LaspyException: Laszip was not found on the system
On some platforms (e.g. Ubuntu 16.04) the installation of laspy does not include laszip/laszip-cli.
See the [issue report](https://github.com/laspy/laspy/issues/79) on github for more infos.

In this case, please follow these steps:

* `wget http://lastools.org/download/LAStools.zip`
* `unzip LAStools.zip && cd LAStools && make`
* `cp bin/laszip /home/USERNAME/miniconda3/envs/pycrown-env/bin/`

If you encounter this error under Windows, please download LAStools.zip, extract the archive and copy the file "laszip.exe" from the "bin"-directory to the conda environment, e.g. C:\Users\<username>\AppData\Local\Continuum\miniconda3\envs\pycrown-env\Scripts\ or C:\Users\<username>\Miniconda3\envs\pycrown-env\Scripts

## Error while building 'pycrown._crown_dalponte_cython' extension
Building the Cython module requires C++ build tools which may need to be installed on your system.

The Windows error message on Windows provides instructions:
`error: Microsoft Visual C++ 14.0 is required. Get it with "Build Tools for Visual Studio": https://visualstudio.microsoft.com/downloads/`
During the setup process, please select 'C++ Build Tools'.

## TypeError: a bytes-like object is required, not 'FakeMmap' when trying to load .laz files
There seems to be an incompatibility between laspy and numpy in recent versions. The combination `numpy==1.16.4` and `laspy==1.5.1` works for me.
I suggest either not using .laz files for the time being or downgrading to the appropiate package versions.
Please also refer to this github issue: https://github.com/laspy/laspy/issues/112


# Getting Started
You can find an IPython Notebook demonstrating each step of the tree segmentation approach in the *example* folder.

You can also run the example python script directly. Results are stored in the *example/result* folder.

`cd example`

`python example.py`

## Main processing steps
### Step 1: Smoothing of CHM using a median filter
![Step 1](example/step_1.jpg)

### Step 2: Tree top detection using local maxima filter
![Step 2](example/step_2.jpg)

### Step 3: Tree Crown Delineation using an adapted Dalponte and Coomes (2016) scheme
![Step 3](example/step_3.jpg)

### Step 4: Tree top correction of trees on steep slopes
![Step 4](example/step_4.jpg)

### Step 5: Smoothing of crown polygons using first returns of normalized LiDAR point clouds
![Step 5](example/step_5.jpg)

### Step 6: Classification of individual trees in the 3D point cloud (visualized with CloudCompare)
![Classified Point Cloud](example/step_6a.jpg)
![Classified Point Cloud](example/step_6b.jpg)
