# Readme.md

This file (README.md) describes the data publication "Data Publication: Refinements in Phase Fraction Determination of Textured Alloys from Transmission Diffraction Data". While readable as a plain text, this file uses Markdown formatting.

Part of the reason to share and document this dataset is to facilitate reuse and expansion, as well as making research more reproducible. Some ideas for expansion include investigating alternate sampling schemes, checking the robustness of the technique with other textures, optimization of sampling schemes, and/or estimation of bias error for a known texture and sampling method.

Questions regarding this data set may be directed to Adam Creuziger (adam.creuziger@nist.gov)

## Sections
- How to Cite
- NIST License
- NIST Disclaimers
- File and Directory Descriptions
- Programs Used
- File Format Descriptions
- Miscellaneous Comments

## Creation Date and Key Updates
Date | Author | Description
--- | --- | ---
2018 Apr 21 | Adam Creuziger | Posted on GitHub
2020 July | Adam Creuziger | Moved JupyterNotebooks to modules, updated to python 3.7 from 2.7
2020 Dec | Adam Creuziger | Collected as companion data to "Refinements in Phase Fraction Determination of Textured Alloys from Transmission Diffraction Data"

# How to Cite

## The data set contents should be cited as:



## This data set is a companion to the following technical paper:



The technical paper describes the motivation and background for this data set, as well as analysis and key insights that were made from this data set.

## Related work

This work builds upon the following papers and data:

Abbreviation | Description
--- | ---
JAC2018 | A. Creuziger, C. Calhoun, W. Poling, T. Gnaupel-Herold "Assessment of Bias Errors Caused by Texture and Sampling Methods in Diffraction-Based Steel Phase Measurements" (2018) Journal of Applied Crystallography, 51, p720-731. DOI: https://doi.org/10.1107/S160057671800420X
JAC2018-data | A. Creuziger, C. Calhoun, W. Poling, T. Gnaupel-Herold, "Data Set: Assessment of Bias Errors Caused by Texture and Sampling Methods in Diffraction-Based Steel Phase Measurements", https://github.com/usnistgov/Texture-Sampling-PhaseMeasurement-BiasErrors
2019Pagan-data |D. Pagan, "Micromechanical Response Quantification using High-energy X-rays during Phase Transformations in Additively Manufactured 17-4 Stainless Steel" Mendeley Data  DOI: 10.17632/3mddz99wsr.1 https://data.mendeley.com/datasets/3mddz99wsr/1


## Citation Guidance

If you are utilizing the results or insights of this work, please cite the technical paper.  If you are utilizing parts of the programs or other documentation provided in the data set, please cite the data set.

# NIST License

See LICENSE.TXT

# NIST Disclaimers

Certain commercial equipment, instruments, or materials are identified in this dataset in order to specify the experimental and computational procedure adequately. Such identification is not intended to imply recommendation or endorsement by the National Institute of Standards and Technology, nor is it intended to imply that the materials or equipment identified are necessarily the best available for the purpose.

# File and Directory Descriptions

Jupyter notebooks, Matlab scripts, and the Rietveld analysis package MAUD constitute the major programs used in this data set. The scripts used by these packages are sorted by directory. When practical, data files generated by these packages are also grouped in the same directory as software that uses them.

This data set is composed of the following:

File or Directory | Description
--- | ---
README.md | File you are currently reading, describing the general structure of the data set.
Workflow.pdf | A flowchart depicting the relationship between figures or tables in the technical paper, input and output files, and scripts or programs.
JupyterNotebooks/ | Directory containing Jupyter notebooks, and associated input and output files.
Matlab/ | Directory containing Matlab scripts, and associated input and output files.
MAUD/ | Directory containing MAUD parameter files, and associated input and output files.

## Workflow.pdf diagram notes (entity-relationship model)

- In the Workflow.pdf diagram, the scripts listed are generally stored in the directory of the program that calls them. The .svg files were created in the program *Inkscape*, which does not have a separate directory.  These .svg files are stored with the image files used to generate the composite figure.

- As the workflow indicates, the Matlab scripts are run first.  Then MAUD is used to create .xpc files.  Then the Jupyter notebooks are run.  The script "PoleFigurePhaseFractions.ipynb" must be run before "PlotTextureCombos.ipynb"; as the averaged intensities must be calculated before they can be plotted.

- To try and minimize crossing lines in the workflow, the script "MaudAnalysis*.par" appears in two locations in the Workflow. The upper appearance is grayed out to indicate the duplication.

- The Matlab script "CalculatedTexturesHWrange.m" does not appear on the Workflow, as this script is nearly identical to "CalculatedTexturesPlots.m".  The "CalculatedTexturesHWrange.m" script saves the plots for a wider range of halfwidths than "CalculatedTexturesPlots.m". Over 150 plots will be generated when running this script.

- Characters '*' and '?' are used to describe folder and file names using the UNIX conventions, where '*' is used for zero or more characters, and '?' is used for exactly as many characters as there are '?' marks.

# Programs Used

### The following application programs and versions were used in this data set:
- Matlab R2015a [1]
 - Mtex 4.3.2 [2]-- with fixes https://github.com/Mtex-toolbox/Mtex/issues/237
- MAUD Version 2.69 (1.0) [3], [4]
- Anaconda 1.3.1, 1.6.2 [5]
- Jupyter Notebook 4.2.3 [6]
- Python 2.7.12
- PF (Peak Fit) [7]
- Inkscape 0.91

The Matlab scripts assume that you have successfully installed the Matlab package Mtex.

### Citations for Programs:
- [1] MATLAB. Natick, MA: MathWorks Inc [Online]. Available: https://www.mathworks.com/product/Matlab.html
- [2] R. Hielscher and H. Schaeben, “A novel pole figure inversion method: specification of the Mtex algorithm,” Journal of Applied Crystallography, vol. 41, no. 6, pp. 1024–1037, Dec. 2008.
- [3] L. Lutterotti and P. Scardi, “Simultaneous structure and size–strain refinement by the Rietveld method,” Journal of Applied Crystallography, vol. 23, no. 4, pp. 246–252, 1990.
- [4] L. Lutterotti, “MAUD - Materials Analysis Using Diffraction.”  [Online]. Available: http://maud.radiographema.eu. [Accessed: 29-Apr-2017]
- [5] Anaconda Software Distribution. Continuum Analytics, 2016 [Online]. Available: https://continuum.io
- [6] F. Perez and B. E. Granger, “IPython: A System for Interactive Scientific Computing,” Computing in Science Engineering, vol. 9, no. 3, pp. 21–29, May 2007.
- [7] T. Gnäupel-Herold, PF. NIST NCNR [Online]. Available: https://www.ncnr.nist.gov/instruments/bt8/BT8DataAnalysis.htm

### Additional Python packages

The environment snapshot (AusteniteModeling.yml) and requirements files (stable-req.txt) listing the python packages used in the jupyter notebook scripts are included in the JupyterNotebooks/ directory. Citations for these additional packages are listed below. With the exception of the mplstereonet package [10], used for plotting pole figures, most of the packages are fairly conventional in scientific computing.

### Citations for additional packages used in jupyter notebooks:

- [8] E. Jones, T. Oliphant, P. Peterson, and others, SciPy: Open source scientific tools for Python. 2001 [Online]. Available: http://www.scipy.org/
- [9] B. Arnold, fortranformat - FORTRAN format interpreter for Python. 2014 [Online]. Available: https://bitbucket.org/brendanarnold/py-fortranformat
- [10] J. Kington, mplstereonet - Stereonets for matplotlib.  [Online]. Available: https://pypi.python.org/pypi/mplstereonet
- [11] fnmatch – Compare filenames against Unix-style glob patterns.  [Online]. Available: https://github.com/python/cpython/blob/2.7/Lib/fnmatch.py
- [12] GitHub.  [Online]. Available: https://github.com
- [13] Glob - Filename globbing utility.  [Online]. Available: https://hg.python.org/cpython/file/2.7/Lib/glob.py
- [14] math - Mathematical functions.  [Online]. Available: https://docs.python.org/2/library/math.html
- [15] os — Miscellaneous operating system interfaces.  [Online]. Available: https://docs.python.org/2/library/os.html
- [16] W. McKinney, “Data Structures for Statistical Computing in Python,” in Proceedings of the 9th Python in Science Conference, 2010, pp. 51–56.
- [17] J. D. Hunter, “Matplotlib: A 2D graphics environment,” Computing In Science & Engineering, vol. 9, no. 3, pp. 90–95, 2007.
- [18] S. van der Walt, S. C. Colbert, and G. Varoquaux, “The NumPy Array: A Structure for Efficient Numerical Computation,” Computing in Science Engineering, vol. 13, no. 2, pp. 22–30, Mar. 2011.

# File Format Descriptions

 File Extension	(Description) | how it was created
 --- | ---
.xlsx (Microsoft Excel Spreadsheet) | generated by a Jupyter Notebook
.pdf (Portable Document Format) |generated by a Jupyter Notebook
.ipynb	(Jupyter/Ipython Notebook)| written using a Jupyter Notebook
.txt	(Plain text)| generated by a Jupyter Notebook, or generated with plain text editor (Atom, Textedit)
.m	(Matlab Code)|  written using Matlab
.png	(Portable Network Graphics)| generated by Matlab or Inkscape
.maa	(BearTex ODF format?)| generated by Matlab.  No explicit format definition for this type of file could be found, so these files were reverse engineered from MAUD output. See notes below in miscellaneous comments.
.sum	(PF pole figure sum file)| generated by PF. This file contains a single line header, and then four columns of whitespace delimited data.
.par	(MAUD analysis file)| generated by MAUD
.prn	(space delimited text file)| generated by Microsoft Excel
.xpc	(BearTex PF format)| generated by MAUD


The programs and files should be platform independent. However, this claim has not been verified.  The dataset was created on Mac OS 10.11.6 and 10.12.5. As such the line endings likely use the line feed (LF) rather than combined carriage return and line feed (CR+LF).

# Miscellaneous Comments

- MAUD uses an .xpc format for pole figures, and .maa format for orientation distribution functions (ODFs), likely derived from BearTex [19]. No explicit documentation of these formats was found. However, upon inspection, this format is similar to the General Intensity File Format in PopLA [20, appendix B2], with a slightly different header and numeric type.  See the following for details:
  - [19] http://eps.berkeley.edu/~wenk/TexturePage/beartex.htm
  - [20] PopLA Manual http://pajarito.materials.cmu.edu/rollett/27750/popLA_Manual.pdf


- The .xpc format reader in PoleFigurePhaseFractions.ipynb was adapted the .epf pole figure reader developed by Youngung Jeong.  The .epf reader is available from:
 - https://github.com/usnistgov/texture/blob/9c0ac8531276a5d8d27c0e895074ca66fda76608/src/upf.py


- As Mtex can also output pole figures, the process of making ODFs in Mtex, exporting the ODFs from Mtex as a .maa file, importing the .maa file into MAUD, and then exporting pole figures in the .xpc format with MAUD may seem redundant and overly complex.  However, the ODFs generated in MAUD were not identical to those created in Mtex when pole figures from Mtex were imported into MAUD. However, when the ODFs were used to generate the textures, the texture indexes matched to the resolution limit of the ODF export. The pole figures also had identical values to the resolution limit of the pole figure export. This was primarily evaluated by comparing the texture index values, as MAUD and Mtex cannot currently plot the ODFs on the same cross section.

- When importing the .maa files into MAUD, select the E-WIMV texture model and 5 degree ODF resolution.  The menu used for import contains a button to export PFs.

- The Matlab script "ExperimentalPFs.m" analyzes and plots orientation distribution functions (ODFs) calculated by pole figure inversion.  The pole figure files are located in the "ExperimentalData/" directory. These pole figures are in a human readable .sum format described above. ODF figures generated by this script are saved in the "ExpFigures/" directory.

- The Matlab script "Calculated_textures.m" creates a series of ODFs with a range of common texture components and a range of texture sharpness.  ODFs in a .maa format readable by the texture analysis program MAUD were saved in the directory "MtexData/", while the phi2=45 (degrees) cross section for each ODF was plotted and saved in the "ODFFigures/" directory.

-  The Matlab scripts "Colormap.m", "ColorMap20.m", "ColorMap4.m", and "ColorMap8.m" were created to plot all ODFs on the same color scale. This color map was created to clearly partition the ODF into values greater than 1 and less than 1, corresponding to the value of a uniform (random) texture distribution.  A grayscale map is used for multiples of a uniform (random) distribution less than 1 and color maps are used for multiples of a uniform (random) distribution equal to or greater than 1. This offers a clear demarcation between regions that have a higher texture volume fraction than uniform and regions that have a lower texture volume fraction than uniform. Using these colormaps require explicitly setting the colormap limits (clim) to match the range used.  For example, when using "ColorMap8.m", the range must be CLim(gcm,[0, 8]), where gcm is the current plot.

- There are many possible settings in MAUD.  Full documentation describing each of these settings in MAUD is under development. When possible, default settings were used.  The .par file contains the parameters used in this analysis and is both human readable and possible to import into MAUD.

- The types of calculations you can perform in MAUD will be extremely limited without a spectra datafile. Files in a generic data format (.prn) can be used. To facilitate calculations, the file "Zeroed-180deg.prn" was created, which is simply a two-column space-delimited list of 2-theta values and intensity values.  In this file the 2-theta values range from 0 to 179.9 degrees in .05 degree increments, and the intensity is set to zero everywhere.


- An internal naming method has been developed by the authors to track individual batches of materials (batch ID) in our research group. The batch IDs listed below uniquely identify a material that comes from the same coil/batch. These batch IDs may be referred to in subsequent publications to define if the material described is the same or a different material.
 - TRIP 780: B080401-AAC-001
 - TRIP 700: B050812-AAC-001

- The simulated intensity profiles (Tilt0Rot0.txt) generated by MAUD for the ND direction (tilt=0 degrees, rotation=0 degrees) for three of the texture combinations (Uniform, TRIP780, TRIP700) are stored in separate directories inside the JupyterNotebooks/ directory. These files were created after loading the texture, simulating the diffraction spectra, plotting, exporting the computed data, and saving as a .cif file.  The header information was removed, data files were renamed and data moved to the JupyterNotebooks/ directory. There is a way to save the files as directly in this format in MAUD, but the authors have been unable to recall or rediscover the steps used to do so.

-  Several of the Jupyter notebooks have additional cells that include information beyond what was included in the paper. These are generally alternative ways to display the data that the authors ultimately chose not to include in the paper. They are included as they may be of use to other researchers.


## Draft-RingRotate Updates

- Additional plots, including difference pole figures, are included in the CHESS_data.m code.
