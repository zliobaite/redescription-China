# Redescription Mining for Ecology and Biogeography, Support material

### Original article:
*Redescription mining for analyzing local limiting conditions: A case study on the biogeography of large mammals in China and southern Asia*. by Esther Galbrun, Hui Tang, Anu Kaakinen and Indrė Žliobaitė, Ecological Informatics, 2021. https://doi.org/10.1016/j.ecoinf.2021.101314

#### Highlights

- We present a methodology for biogeographical analysis
- Redescription mining emphasizes local association patterns and limiting conditions
- Redescription mining combines different perspectives over the studied system
- We showcase the potential of this method for ecological and biogeographical studies
- We consider an example biogeographic study focused on China and southern Asia

#### Abstract

Identifying and understanding limiting conditions is at the centre of ecology and biogeography.
Traditionally, associations between climate and occurrences of organisms are inferred from observational data using regression analysis, correlation analysis or clustering.
Those methods extract patterns and relationships that hold throughout a dataset.
We present a computational methodology called redescription mining, that emphasizes local patterns and associations that hold strongly on subsets of the dataset, instead.
We aim to showcase the potential of this methodology for ecological and biogeographical studies, and encourage researchers to try it.

Redescription mining can be used to identify associations between different descriptive views of the same system. It produces an ensemble of local models, that provide different perspectives over the system.
Each model (redescription) consists of two sets of limiting conditions, over two different views, that hold locally. Limiting conditions, as well as the corresponding subregions, are identified automatically using data analysis algorithms.

We explain how this methodology applies to a biogeographic case study focused on China and southern Asia.
We consider dental traits of the large herbivorous mammals that occur there and climatic conditions as two aspects of this ecological system, and look for associations between them.

Redescription mining can offer more refined inferences on the potential relation between variables describing different aspects of a system than classical methods. Thus, it permits different questions to be posed of the data, and can usefully complement classical methods in ecology and biogeography to uncover novel biogeographic patterns.

A python package for carrying out redescription mining analysis is publicly available.

## List of contents

- The **data** folder contains the datasets used in the study.
- The **python-siren** folder contains a copy of the ReReMi and Siren source code (v4.3.0, python 2) as used to carry out the redescription mining analysis.
- The **standard-methods** contains the scripts used to carry out the analysis with standard methods, i.e. correlation analysis, principal component analysis (PCA), regression analysis and clustering.
- The **figures** folder contains additional high-resolution figures.

## Running the analysis with standard methods

A python notebook containing the script for carrying out the analysis with standard methods, is provided in folder *standard-methods* (`methods_compare.ipynb`). The output of the notebook has been exported to html to allow viewing the results without having to re-run the analysis (`methods_compare.html`). 


## Running the analysis with redescription mining

1. **Obtain the tools for redescription mining.** 
    Use either the code provided in folder *python-siren* or the latest version, obtained following instructions from the [website](http://cs.uef.fi/siren/main/download.html). Note that some packages need to be installed for the code to work.

2. **Mine the datasets.**

    
    First, create a folder for the output

    ```
    mkdir ./xps
    ```

    launch the first run

    ```
    python PATH_TO_PYTHON-SIREN/siren/reremi/mainReReMi.py preferences_i.01o.3CONJ2.xml
    ```

    and the second

    ```
    python PATH_TO_PYTHON-SIREN/siren/reremi/mainReReMi.py preferences_i.01o.3m.1.xml
    ```

    This will generate files containing redescriptions (`.queries`) and log (`.txt`) in the `xps` folder.

3. **Load the redescription in Siren.**

    It is possible to load a list of redescriptions into Siren with the following command

    ```
    python PATH_TO_PYTHON-SIREN/exec_siren.py preferences_i.01o.3CONJ2.xml
    ```

    and correspondingly for the other run.

    Alternatively, a package containing the full lists of redescriptions and the shortlist presented in the study, as well as adjusted parameters for filtering and plotting is provided and can be loaded into Siren with the following command

    ```
    python PATH_TO_PYTHON-SIREN/exec_siren.py biotraits_focus-sel.siren
    ```

4. **Filter redescriptions.**
    The obtained lists of redescriptions where filtered separately to remove redundant redescriptions. This is done by setting the `max_overlap_rows` parameter (`Menu File > Preferences... > Filtering`) to 0.9, then going to the redescription tab, selecting the chosen list of redescription, ensuring that it is sorted by decreasing redescription accuracy (clicking in the `J` column) then right-clicking on the first redescription and selecting `Filter row redundant pairwise` from the pop-up menu.

5. **Visualize redescriptions.**
    Individual redescriptions can be visualized in various ways (including maps and parallel coordinates plots) by right-clicking on the redescription of interest, and selecting `View` and the desired type of visualization.

6. **Visualize redescription summaries.**
    A visual summary of a list of redescriptions, produced by clustering entities depending on which redescriptions they support, can be generated by right-clicking on the redescription list of interest, and selecting `View` and the desired type of visualization.
