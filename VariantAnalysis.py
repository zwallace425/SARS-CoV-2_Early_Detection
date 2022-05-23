import re
import pandas as pd
import numpy as np


# Compute the Sequence Prevalence Score for single amino acid substitution either within a specific domain,
# the entire spike protein, or for a user inputted list.  When ranking substitutions within a user inputted list,
# the algorithm will simply rank the substitutions and not consider the domain.
def mutation_ranking(variant_df, interval, domain, mutations = []):
	if (mutations):
		variant_df = variant_df[(variant_df['Variant'].isin(mutations))]
	elif (domain == "NTD"):
		variant_df = variant_df[(variant_df['Position'] >= 13) & (variant_df['Position'] <= 303)]
	elif (domain == "RBD"):
		variant_df = variant_df[(variant_df['Position'] >= 319) & (variant_df['Position'] <= 541)]
	elif (domain == "Other"):
		variant_df = variant_df[(variant_df['Position'] >= 542) & (variant_df['Position'] <= 1273)]	

	variant_df = variant_df[variant_df[variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,0].name] > 10]
	recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:interval]
	recent_growth_rates = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:interval]
	mut_growth_data = []
	mut_prevalence_data = []
	for i in range(interval):
		growth_data = pd.DataFrame({'Variant': variant_df['Variant'], 'Growth Rate': recent_growth_rates[recent_growth_rates.iloc[:,i].name]})
		prevalence_data = pd.DataFrame({'Variant': variant_df['Variant'], 'Prevalence': recent_prevalence[recent_prevalence.iloc[:,i].name]})
		mut_growth_data.append(growth_data)
		mut_prevalence_data.append(prevalence_data)
	mut_growth_data = pd.concat(mut_growth_data, ignore_index = True)
	mut_prevalence_data = pd.concat(mut_prevalence_data, ignore_index = True)
	mut_data = pd.DataFrame(data = {'Variant': mut_growth_data['Variant'], 'Growth Rate': mut_growth_data['Growth Rate'], 'Prevalence': mut_prevalence_data['Prevalence']})
	mut_data = mut_data.dropna()
	significant_mutants = mut_data[(mut_data['Growth Rate'] > 5) | (mut_data['Prevalence'] > 0.05)]
	mutation_ranking = significant_mutants.groupby('Variant').size().reset_index(name = 'Mutation Prevalence Score').sort_values(by = 'Mutation Prevalence Score', ascending = False)
	prevalence_median = significant_mutants.groupby(['Variant'])['Prevalence'].median().reset_index(name = 'Prevalence Median')
	mutation_ranking = mutation_ranking.merge(prevalence_median, on = "Variant")
	mutation_ranking = mutation_ranking.sort_values(by = ['Mutation Prevalence Score','Prevalence Median'], ascending = [False, False]).drop_duplicates().reset_index(drop = True)
	return(mutation_ranking)

# Compute the Composite Score by summing the Substitution Ranking score with the Functional Ranking score and then
# ranking the results to get the 'Overall Spike Rank'.  This is then then a combined score for interpreting
# the overall threat level of an emerging variant.
def composite_ranking(variant_df, sfoc_df, interval, who = [], lineage = [], covariants = []):
	sequence_prevalence_score = sequence_ranking(variant_df, interval, who, lineage, covariants)
	functional_impact_score = functional_ranking(variant_df = pd.DataFrame({'Variant': sequence_prevalence_score['Variant']}), sfoc_df = sfoc_df, covariants = list(sequence_prevalence_score['Variant']))
	composite_score = pd.merge(sequence_prevalence_score, functional_impact_score, on = "Variant")
	composite_score['Composite Score'] = composite_score['Sequence Prevalence Score'] + composite_score['Functional Impact Score']
	overall_spike_rank = composite_score.sort_values(by = ['Composite Score', 'Prevalence Median'], ascending = [False, False]).drop_duplicates().reset_index(drop=True)
	return(overall_spike_rank)


# Compute the Sequence Prevalence Scores for a specific lineage as described on slide 15 of
# Richard's early detection presintation.  This basically ranks spike covariants by considering the
# sequence prevalence dynamics, growth rates, and geographic spread by region.  This algorithm essentially 
# ranks individual covariants based on their epedimological dynamics or lineage expansions.
def sequence_ranking(variant_df, interval, who = [], lineage = [], covariants = []):
	all_data = True
	if (who):
		variant_df = variant_df[(variant_df['WHO Label'].isin(who))]
		all_data = False
	if (lineage):
		variant_df = variant_df[(variant_df['PANGO Lineage'].isin(lineage))]
		all_data = False
	if (covariants):
		variant_df = variant_df[(variant_df['Variant'].isin(covariants))]
		all_data = False	
	variant_df = variant_df[variant_df[variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,0].name] > 10]
	recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:interval]
	recent_growth_rates = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:interval]
	variant_growth_data = []
	variant_prevalence_data = []
	for i in range(interval):
		growth_data = pd.DataFrame({'WHO Label': variant_df['WHO Label'], 'Variant': variant_df['Variant'], 'Growth Rate': recent_growth_rates[recent_growth_rates.iloc[:,i].name]})
		prevalence_data = pd.DataFrame({'WHO Label': variant_df['WHO Label'], 'Variant': variant_df['Variant'], 'Prevalence': recent_prevalence[recent_prevalence.iloc[:,i].name]})
		variant_growth_data.append(growth_data)
		variant_prevalence_data.append(prevalence_data)
	variant_growth_data = pd.concat(variant_growth_data, ignore_index = True)
	variant_prevalence_data = pd.concat(variant_prevalence_data, ignore_index = True)
	variant_data = pd.DataFrame(data = {'WHO Label': variant_growth_data['WHO Label'], 'Variant': variant_growth_data['Variant'], 'Growth Rate': variant_growth_data['Growth Rate'], 'Prevalence': variant_prevalence_data['Prevalence']})
	variant_data = variant_data.dropna()
	significant_variants = variant_data[(variant_data['Growth Rate'] > 5) | (variant_data['Prevalence'] > 0.05)]
	who_and_variants = significant_variants[["WHO Label", "Variant"]].drop_duplicates()
	covariant_ranking = significant_variants.groupby(['Variant']).size().reset_index(name = 'Sequence Prevalence Score')
	prevalence_median = significant_variants.groupby(['Variant'])['Prevalence'].median().reset_index(name = 'Prevalence Median')
	covariant_ranking = covariant_ranking.merge(prevalence_median, on = "Variant")
	unscored_variants = sorted(set(covariants) - set(list(covariant_ranking['Variant'])), key = covariants.index)
	unscored_variant_df = pd.DataFrame({'Variant': unscored_variants, 'Sequence Prevalence Score': [0]*len(unscored_variants), 'Prevalence Median': [0]*len(unscored_variants)})
	covariant_ranking = pd.concat([covariant_ranking, unscored_variant_df], axis = 0)
	if (all_data):
		covariant_ranking = who_and_variants.merge(covariant_ranking, on = 'Variant')
	covariant_ranking = covariant_ranking.sort_values(by = ['Sequence Prevalence Score','Prevalence Median'], ascending = [False, False]).drop_duplicates().reset_index(drop = True).head(50)
	return(covariant_ranking)


# Compute the Functional Score for a covariant based on the defined criteria for an SFoC.
# Return a dataframe with the covariant and the Predcited Functional Score, ranked by the score.
# This takes in the "World - Variants" tab extracted from the Emerging Variants Report and a the
# "SFoCs" tab extracted from the Emerging Variants Report.  The SFoCs tab is a pre-defined table
# of "Sequence Features of Concern", where the JCVI team manualy curated a list of sequence regions
# based on experimental data to be defined as regions that could have a functional consequence on 
# the spike, such antibody neutralization.
def functional_ranking(variant_df, sfoc_df, who = [], lineage = [], covariants = []):
	if (who):
		variant_df = variant_df[(variant_df['WHO Label'].isin(who))]
	if (lineage):
		variant_df = variant_df[(variant_df['PANGO Lineage'].isin(lineage))]
	variant_df.drop_duplicates(subset = ['Variant']).reset_index(inplace = True, drop = True)
	covariants = list(variant_df['Variant'])
	scores = []
	print('Computing Functional Impact Scores ...')
	for i in range(len(covariants)):
		mutations = covariants[i].split(",")
		sfoc_score = 0
		prev_pos = 1
		prev_mut = "N"
		for j in range(len(mutations)):
			aa_pos = int(re.sub("[^0-9]", "", mutations[j]))
			aa_muts = re.sub(r'[0-9]+', '', mutations[j])
			mut = aa_muts[len(aa_muts) - 1]
			sfocs = sfoc_df[(sfoc_df['Start'] <= aa_pos) & (sfoc_df['End'] >= aa_pos)]
			sfocs.reset_index(inplace = True, drop = True)
			if (not ((aa_pos == (prev_pos + 1)) and (mut == "-") and (prev_mut	== "-"))):
				for k in range(len(sfocs)):
					if (aa_pos == 614):
						sfoc_score += 1
					if (pd.notna(sfocs.at[k,"mAb escape"])):
						if ("class 1" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 2" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 3" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
						if ("class 4" in sfocs.loc[k].at["mAb escape"]):
							sfoc_score += 1
					if (pd.notna(sfocs.at[k,"serum Ab escape"])):
						if ("convalescent serum" in sfocs.loc[k].at["serum Ab escape"]):
							sfoc_score += 1
						if ("Moderna vaccine serum" in sfocs.loc[k].at["serum Ab escape"]):
							sfoc_score += 1
					if (pd.notna(sfocs.at[k,"Increased ACE2 binding"])):
						sfoc_score += 1
					if (pd.notna(sfocs.at[k,"Region of interest"])):
						sfoc_score += 1
			prev_pos = aa_pos
			prev_mut = mut
		scores.append(sfoc_score)
	
	functional_impact = pd.DataFrame({'Variant': covariants, 'Functional Impact Score': scores})
	impact_ranking = functional_impact.sort_values(by = 'Functional Impact Score', ascending = False)
	impact_ranking = impact_ranking.reset_index(drop=True)
	return(impact_ranking)


# Compute the Emergance Score for a PANGO lineage using the Emerging Lineage Ranking heuristic.
# Return a dataframe with the lineage and emergence score, ranked by the Emergence Score.
# This takes in the "World - Variants" tab extracted from the Emerging Variants Report and
# returns a score for lineages with at least 10 variants counts in the recent months and meets
# a significant growth rate level.  Note, this heuristic is not used very  much in practice.
def lineage_ranking(variant_df, interval):
	variant_df = variant_df[variant_df[variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,0].name] > 10]
	recent_growth_rates = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:interval]
	recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:interval]
	who_and_pango = variant_df[["WHO Label", "PANGO Lineage"]].drop_duplicates()
	lineage_growth_data = []
	lineage_prevalence_data = []
	for i in range(interval):
		growth_data = pd.DataFrame({'WHO Label': variant_df['WHO Label'], 'PANGO Lineage': variant_df['PANGO Lineage'], 'Growth Rates': recent_growth_rates[recent_growth_rates.iloc[:,i].name]})
		prevalence_data = pd.DataFrame({'WHO Label': variant_df['WHO Label'], 'PANGO Lineage': variant_df['PANGO Lineage'], 'Prevalence': recent_prevalence[recent_prevalence.iloc[:,i].name]})
		lineage_growth_data.append(growth_data)
		lineage_prevalence_data.append(prevalence_data)
	lineage_growth_data = pd.concat(lineage_growth_data, ignore_index = True)
	lineage_prevalence_data = pd.concat(lineage_prevalence_data, ignore_index = True)
	lineage_data = pd.DataFrame(data = {'WHO Label': lineage_growth_data['WHO Label'], 'PANGO Lineage': lineage_growth_data['PANGO Lineage'], 'Growth Rates': lineage_growth_data['Growth Rates'], 'Prevalence': lineage_prevalence_data['Prevalence']})
	lineage_data = lineage_data[(lineage_data['Growth Rates']) > 5 | (lineage_data['Prevalence'] > 0.05)]
	emergence_ranking = lineage_data.groupby('PANGO Lineage').size().reset_index(name = 'Emergence Score')
	emergence_ranking = who_and_pango.merge(emergence_ranking, on = 'PANGO Lineage').sort_values(by = 'Emergence Score', ascending = False).reset_index(drop=True)
	return(emergence_ranking)