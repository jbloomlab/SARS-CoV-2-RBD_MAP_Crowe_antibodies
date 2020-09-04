"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os.path
import textwrap
import urllib.request

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs.groupby(['library', 'sample'])) == len(barcode_runs)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *selection*.
sample_vs_expect = (
    barcode_runs
    .assign(expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'selection']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect['equal'].all(), sample_vs_expect.query('equal != True')

# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        codon_variant_table=config['codon_variant_table'],
        count_variants=nb_markdown('count_variants.ipynb'),
        variant_counts=config['variant_counts'],
        counts_to_scores=nb_markdown('counts_to_scores.ipynb'),
        scores_to_frac_escape=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
        escape_fracs_homologs=config['escape_fracs_homologs'],
        analyze_escape_profiles=nb_markdown('analyze_escape_profiles.ipynb'),
        mds_escape_profiles=nb_markdown('mds_escape_profiles.ipynb'),
        significant_escape_sites=config['significant_escape_sites'],
        output_pdbs='results/summary/output_pdbs.md',
        circulating_variants='results/summary/circulating_variants.md',
        crowe_evol='results/summary/evolution_escape_Crowe.md',
        make_supp_data=nb_markdown('make_supp_data.ipynb'),
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each notebook in the workflow:

            1. Get codon-variant-table from [here]({config['codon_variant_table_url']}).

            2. [Count variants]({path(input.count_variants)}) to create a
               [variant counts file]({path(input.variant_counts)}).

            3. [Escape scores from variant counts]({path(input.counts_to_scores)}).

            4. [Escape fractions for mutations and homologs]({path(input.scores_to_frac_escape)});
               creating [mutation escape fraction file]({path(input.escape_fracs)})
               and [homolog escape fraction file]({path(input.escape_fracs_homologs)}).

            5. Analyze and plot [escape profiles]({path(input.analyze_escape_profiles)}).
               Also write file with [sites of strong escape]({path(input.significant_escape_sites)}).

            6. [Multidimensional scaling]({path(input.mds_escape_profiles)}) on escape profiles.

            7. Map escape profiles to ``*.pdb`` files using [this notebook]({path(input.output_pdbs)})

            8. [Create list of circulating variants from GISIAID sequences]({path(input.circulating_variants)}).

            9. Analyze escape profiles in light of evolutionary and functional constraint using [this notebook]({path(input.crowe_evol)}).

            10. [Make supplementary data files]({path(input.make_supp_data)}),
                which are [here]({path(config['supp_data_dir'])}).

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule make_supp_data:
    input:
        config['escape_profiles_config'],
        config['escape_fracs'],
        config['escape_profiles_dms_colors']
    output:
        nb_markdown=nb_markdown('make_supp_data.ipynb'),
        outdir=directory(config['supp_data_dir']),
    params:
        nb='make_supp_data.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule crowe_evol:
    input:
        config['significant_escape_sites'],
        config['circulating_variants']
    output:
        md='results/summary/evolution_escape_Crowe.md',
        md_files = directory('results/summary/evolution_escape_Crowe_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='evolution_escape_Crowe.Rmd',
        md='evolution_escape_Crowe.md',
        md_files='evolution_escape_Crowe_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule output_pdbs:
    input:
        config['escape_fracs']
    output:
        md='results/summary/output_pdbs.md'
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='output_pdbs.Rmd',
        md='output_pdbs.md'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md}
        """

rule circulating_variants:
    input:
        config['significant_escape_sites']
    output:
        md='results/summary/circulating_variants.md',
        table = config['circulating_variants']
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='circulating_variants.Rmd',
        md='circulating_variants.md'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md}
        """

rule mds_escape_profiles:
    """Multi-dimensional scaling on antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs'],
        escape_profiles_config=config['escape_profiles_config'],
        site_color_schemes=config['site_color_schemes'],
    output:
        nb_markdown=nb_markdown('mds_escape_profiles.ipynb'),
    params:
        nb='mds_escape_profiles.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule analyze_escape_profiles:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs'],
        escape_profiles_config=config['escape_profiles_config'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence'],
        mut_bind_expr=config['mut_bind_expr'],
    output:
        nb_markdown=nb_markdown('analyze_escape_profiles.ipynb'),
        significant_escape_sites=config['significant_escape_sites'],
        escape_profiles_dms_colors=config['escape_profiles_dms_colors'],
    params:
        nb='analyze_escape_profiles.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule scores_to_frac_escape:
    """Estimate mutation- and homolog-level escape scores."""
    input:
        escape_score_samples=config['escape_score_samples'],
        escape_scores=config['escape_scores'],
        escape_scores_homologs=config['escape_scores_homologs'],
    output:
        nb_markdown=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
        escape_fracs_homologs=config['escape_fracs_homologs'],
    params:
        nb='scores_to_frac_escape.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule counts_to_scores:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts'],
        config['wildtype_sequence'],
        config['mut_bind_expr'],
        config['variant_expr'],
        config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores.ipynb'),
        escape_scores=config['escape_scores'],
        escape_scores_homologs=config['escape_scores_homologs'],
        escape_score_samples=config['escape_score_samples'],
    params:
        nb='counts_to_scores.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule count_variants:
    """Count variants from Illumina barcode runs."""
    input:
        config['codon_variant_table'],
        config['barcode_runs'],
        config['wildtype_sequence']
    output:
        config['variant_counts'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_codon_variant_table:
    """Download codon variant table from URL."""
    output:
        codon_variant_table=config['codon_variant_table']
    run:
        urllib.request.urlretrieve(config['codon_variant_table_url'],
                                   output.codon_variant_table)

rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)

rule get_variant_expr:
    """Download SARS-CoV-2 variant expression from URL."""
    output:
        file=config['variant_expr']
    run:
        urllib.request.urlretrieve(config['variant_expr_url'], output.file)

rule get_variant_bind:
    """Download SARS-CoV-2 variant ACE2-binding from URL."""
    output:
        file=config['variant_bind']
    run:
        urllib.request.urlretrieve(config['variant_bind_url'], output.file)
