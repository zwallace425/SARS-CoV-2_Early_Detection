README

DEPENDENCIES: python >= 3.x, pandas

TITLE: SARS-CoV-2 Spike Protein Variant Analysis Automation Pipeline

INTRODUCTION: This software packages uses all the data that has been processed through BV-BRC GISAID pipeline and provides extensive analysis.
BV-BRC collects data has epedimiological data on emerging Spike protein variants including lineage counts, prevalence, and growth by country,
covariant counts, prevalence, and growth by country, and single point mutations counts, prevalence, and growth by country.  All
that data is processed from GISAID on a weekly basis.  To offer a range of methods for analyzing growth dynamics of emerging SARS-CoV-2
variant, we this offer this analysis pipeline.  It deploys a range of heuristic methods for ranking variants based on sequence prevalence
dynamics as well as functional impact dynamics.  So basically this pipeline offers ways to prioritize variants with score and also visualize
growth of variants with plots.

FILES: This pipeline requires the Emerging Variants Report file generated from Mauliks pipeline.  There are also other optional files that
can be used in the analysis.

USAGE: Essentially there are two primary usages of this pipeline: ranking and graphing.  The ranking algorithms could either rank PANGO
lineages using the Emerging Lineage Ranking algorithm, rank covariants/constellations of variants using the Substitution Ranking, Functional
Ranking, or Composite Ranking algorithm, or rank single point mutations using the Mutation Ranking Algorithm.

To run the pipeline, try running the following first to see the usage to the commandline:

	python main.py --filename [Emerging Variants Report] --analysis help

This will print out both the ranking usage and the graph usage.  Please follow that thoroughly as those commands are the only way to run
the software. Any mistake and the sofware will terminate and report those usage messages again.

The '--analysis' commandline argument must always be present to run this software. The following are the permitted options for the
'--analysis':

	emergence_ranking, substitution_ranking, functional_ranking, composite_ranking, mutation_ranking, graph, help

There are other madatory arguments following some of these analysis options, but use the 'help' to see this.