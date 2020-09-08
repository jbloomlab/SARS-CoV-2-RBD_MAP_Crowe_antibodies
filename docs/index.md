# High-resolution mapping of the neutralizing and binding specificities of polyclonal rabbit serum elicited by HIV trimer immunizations

Below are links to interactive visualizations of the HIV mutational antigenic profiling of sera from rabbits vaccinated with BG505 trimers, enabled by [dms-view](https://dms-view.github.io/docs/). These data are posted in [this pre-print - broken link](link).

[Figure 2 - Plotting both the glycan hole and C3/V5 epitopes](https://dms-view.github.io/?markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.csv&condition=2124-Wk22&site_metric=site_positive+diffsel&mutation_metric=mut_pos+mutdiffsel&selected_sites=84%2C85%2C86%2C87%2C88%2C89%2C90%2C229%2C230%2C231%2C240%2C241%2C242%2C243%2C268%2C289%2C290%2C291%2C347%2C350%2C351%2C352%2C353%2C354%2C355%2C356%2C357%2C358%2C359%2C360%2C396%2C459%2C460%2C461%2C462%2C463%2C464%2C465%2C466%2C467%2C629&pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2F5fyl_trimer_renumber.pdb)

[Figure 4 - Plotting the C3/V5 epitope](https://dms-view.github.io/?markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.csv&condition=5724-Wk26&site_metric=site_positive+diffsel&mutation_metric=mut_pos+mutdiffsel&selected_sites=350%2C351%2C352%2C353%2C354%2C355%2C356%2C357%2C358%2C359%2C360%2C396%2C459%2C460%2C461%2C462%2C463%2C464%2C465%2C466%2C467&pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2F5fyl_trimer_renumber.pdb)

[Figure 5 - Plotting the glycan hole epitope](https://dms-view.github.io/?markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2FDingens2020.csv&condition=2124-Wk22&site_metric=site_positive+diffsel&mutation_metric=mut_pos+mutdiffsel&selected_sites=84%2C85%2C86%2C87%2C88%2C89%2C90%2C229%2C230%2C231%2C240%2C241%2C242%2C243%2C268%2C289%2C290%2C291%2C347%2C629&pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FHIV%2Fmaster%2Fdata%2FEnv%2FDingens2020%2F5fyl_trimer_renumber.pdb)

For dms-view visualizations, logo plots are colored as in Figure 2 of the paper. The glycan hole epitope is green, the C3/V5 epitope is blue, and all other sites are grey. These residue-level epitope annotations are arbitrarily defined based on all data, we encourage readers to further explore additional sites. Mutations tested in preliminary point mutant mapping (Figure 1A) are shown in black. Selection is mapped onto the [5fyl BG505 trimer or monomer structure](https://www.rcsb.org/structure/5FYL).

[This page](https://jbloomlab.github.io/dms_tools2/diffsel.html) documents the differential selection statsitics we use to analyze these data.

The site-metrics (dot plot) include:

- **positive diffsel**: The sum of all positive differential selection values at a site. This gives a sense to the total amount of escape/selective pressure at each site.
- **negative diffsel**: The sum of all negative differential selection values at a site. This gives a sense for mutations that are depleted, rather than enriched, during serum selection relative to a non-selected control library. It is intriguing that many of these potential serum sensitizing mutations cluster and are consistent across sera.
- **max diffsel**: The value of the largest effect mutation (largest mutation differential selection) at each site. This gives a sense of the maximal effect of any mutation at each site.
- **min diffsel**: The value of the smallest effect mutation (smalled mutation differential selection) at each site.
- **abs diffsel**: The absolute value of all mutation differential selection values at each site.



The mutation-metrics (logoplot) include

- **diffsel**: All mutation differential selection values, including negative values, are plotted.
- **pos diffsel**: Only the mutations with positive differential selection values.

Additionally,

- The frequency at which each amino acid is found in nature (**Natural Frequencies**), accessed from [LANL's filtered web alignment](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html]) (2017 version).
- The BG505 amino-acid preferences (**DMS preferences**), determined using the same BG505.T332N mutant virus libraries in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). Here, the height of each amino acid is proportional to how well that virus replicates in cell culture. This statistic can crudely be used to examine what mutations are viable and in our mutant virus libraries before serum selection.