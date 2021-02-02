# Complete mapping of mutations to the SARS-CoV-2 spike receptor-binding domain that escape antibody recognition
See the [paper describing the work here](https://www.sciencedirect.com/science/article/pii/S1931312820306247).

Study and analysis by Allie Greaney, Tyler Starr, Pavlo Gilchuk, Seth Zost, ..., [James Crowe](https://www.vumc.org/crowe-lab) and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html), and co-authors.

## Summary of workflow and results
For a summary of the workflow and links to key results files, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment.

 2. The required input data.

 3. The computer code and a [Snakemake](https://snakemake.readthedocs.io) file to run it.


### Build the computing environment

First, set up the computing environment, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [environment.yml](environment.yml).
If you have not previously built the conda environment, then build it to `./env` with:

    conda env create -f environment.yml -p ./env

Then activate it with:

    conda activate ./env

If you've previously built the environment into `./env`, just do the activation step.

Setting up the `conda` environment above installs everything to run all parts of the analysis **except** the `R` markdown notebooks.
For those, the pipeline currently uses the Fred Hutch computing cluster module `R/3.6.1-foss-2018b` as specified in `Snakefile`.
That module is not packaged with this repo, so if you aren't on the Fred Hutch cluster you'll have to create a similar `R` environment yourself (all the `R` packages are listed at the beginning of their output in the [summary results](results/summary/summary.md).

### Input data
The input data are specified in [./data/](data); see the README in that subdirectory for more details.

### Running the code
The analysis consists of Jupyter and R Markdown notebooks in the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment in `./env`, as in:

    snakemake --use-conda --conda-prefix ./env --use-envmodules -R make_summary -j 1

However, you probably want to using a cluster to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the Fred Hutch server, simply run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the Hutch server resources.
You likely want to submit [run_Hutch_cluster.bash](run_Hutch_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 7-0 run_Hutch_cluster.bash

## Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml).
This file defines key variables for the analysis, and should be relatively self-explanatory.
You should modify the analysis by changing this configuration file; do **not** hard-code crucial experiment-specific variables within the notebooks or `Snakefile`.

## Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the Fred Hutch cluster.
The [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script uses this configuration to run [Snakefile](Snakefile).
If you are using a different cluster than the Fred Hutch one, you wll need to modify the cluster configuration file.

## Notebooks that perform the analysis
The Jupyter notebooks and R Markdown scripts that perform most of the analysis are in this top-level directory with the extension `*.ipynb` or `*.Rmd`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

## Results
Results are placed in the [./results/](results) subdirectory.
Many of the files created in this subdirectory are not tracked in the `git` repo as they are very large.
However, key results files are tracked as well as a summary that shows the code and results.
Click [here](./results/summary/summary.md) to see that summary.

The large results files are tracked via [git-lfs](https://git-lfs.github.com/).
This requires `git-lfs` to be installed, which it is in the `conda` environment specified by [environment.yml](environment.yml).
The following commands were then run:

    git lfs install

You may need to run this if you are tracking these files and haven't installed `git-lfs` in your user account.
Then the large results files were added for tracking with:

    git lfs track <FILENAME>

## Updating the conda environment
[environment.yml](environment.yml) contains a fully pinned conda environment.
An environment without all of the versions pinned is in [environment_unpinned.yml](environment_unpinned.yml).
If you need to update the environment, the suggested way to do it is add the new requirement to [environment_unpinned.yml](environment_unpinned.yml), then build that environment and finally export the pinned version:

    conda env create -f environment_unpinned.yml --prefix ./env

Then activate the environment with:

    conda activate ./env

Finally, export the pinned version with:

    conda env export --prefix ./env > environment.yml

## Creating "subset" repos and uploading data to the SRA
Currently this repo contains analyses of many antibodies and sera, and should remain public since collaborators do not want all of these data to be public.

For papers, you can make a public "subset" repo by following the instructions in [./subset_data/](subset_data).
After making a subset repo, you can upload sequencing data to the Sequence Read Archive (SRA) following the instructions in [./SRA_upload/](SRA_upload).
