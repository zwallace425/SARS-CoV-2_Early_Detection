import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import seaborn as sns
from functools import reduce
from operator import itemgetter




# This function is meant to plot the prevalence ratios (as stacked plot) or the growth rate trends of single amino acid
# point mutations either in a user specified domain (NTD, RBD, other), the entire spike protein, or a user inputted list.
# For simplicity, this will only plot the 10 mutations with the highest prevalence or growth rates accross the past 6 months
def plot_mutations(variant_df, region, domain, mutations = []):
	graph_type = input("Plot the changes over time of the mutation's prevalence or growth rate in "+region+"? [prevalence/growth]: ")
	if (mutations):
		variant_df = variant_df[(variant_df['Variant'].isin(mutations))]
	if (domain == "NTD"):
		variant_df = variant_df[(variant_df['Position'] >= 13) & (variant_df['Position'] <= 303)]
	elif (domain == "RBD"):
		variant_df = variant_df[(variant_df['Position'] >= 319) & (variant_df['Position'] <= 541)]
	elif (domain == "Other"):
		variant_df = variant_df[(variant_df['Position'] >= 542) & (variant_df['Position'] <= 1273)]	
		domain = "'Other Domain'"

	if (graph_type == 'prevalence'):
		recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:6]
		recent_prevalence = recent_prevalence.fillna(0)
		mutation_prevalence_df = pd.concat([pd.DataFrame({'Variant': variant_df['Variant']}), recent_prevalence], axis = 1)
		mutation_prevalence_df = mutation_prevalence_df.reindex(sorted(mutation_prevalence_df.columns), axis = 1)
		mutation_prevalence_df = mutation_prevalence_df.set_index('Variant')
		mutation_prevalence_df['rowmax'] = mutation_prevalence_df.max(axis = 1)
		mutation_prevalence_df = mutation_prevalence_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		mutation_prevalence_df.columns = mutation_prevalence_df.columns.str.replace("Prevalence - ", "")
		x_val = np.asarray(mutation_prevalence_df.columns.tolist())
		indices = mutation_prevalence_df.index.tolist()
		for i in indices:
			mask = np.isfinite(mutation_prevalence_df.loc[i].values)
			plt.plot(x_val[mask], mutation_prevalence_df.loc[i].values[mask], label = i)
		plt.suptitle('A Time Course of Prevalence Ratios for '+domain+' Mutations from Sequences Isolated in '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Mutation Prevalence Ratio (Mutation Count/Isolates Count)')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()

	elif (graph_type == 'growth'):
		recent_growth_rate = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:6]
		recent_growth_rate = recent_growth_rate.fillna(0)
		mutation_growth_df = pd.concat([pd.DataFrame({'Variant': variant_df['Variant']}), recent_growth_rate], axis = 1)
		mutation_growth_df = mutation_growth_df.reindex(sorted(mutation_growth_df.columns), axis = 1)
		mutation_growth_df = mutation_growth_df.set_index('Variant')
		mutation_growth_df['rowmax'] = mutation_growth_df.max(axis = 1)
		mutation_growth_df = mutation_growth_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		mutation_growth_df.columns = mutation_growth_df.columns.str.replace("Growth Rate - ", "")
		x_val = np.asarray(mutation_growth_df.columns.tolist())
		indices = mutation_growth_df.index.tolist()
		for i in indices:
			mask = np.isfinite(mutation_growth_df.loc[i].values)
			plt.plot(x_val[mask], mutation_growth_df.loc[i].values[mask], label = i)
		plt.suptitle('A Time Course of Growth Rate Trends for '+domain+' Mutations from Sequences Isolated in '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Growth Rate')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()
	else:
		sys.exit("Invalid graph type entry. Please try again and correctly specify the appropriate graph type (prevalence or growth)")


# This is a helper function used for plot_covariants() meant to assist with plotting covariants that are duplicated
# within the variants_df, which arises when the variant_df is unfiltered, filtered for WHO label, or filtered for 
# specific covariants.  The only case when covariants are not duplicated through a dataframe is when it is filtered
# by the pango lineage.  Since in those cases mentioned the covariants we wish to plot are duplicated, the prevalence
# and growth rate dynamics need to be recalculated to reflect the sum of all variant counts for that particular
# covariant.
def plot_covariants_help(variant_df, covariants):
	
	def div(x, y):
		if y == 0:
			return 0
		return x/y	

	prevalence_cols = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:6]
	prevalence_cols = prevalence_cols.columns
	growth_cols = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:6]
	growth_cols = growth_cols.columns
	variant_cols = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,:6]
	variant_cols = variant_cols.columns
	new_data = []
	for cov in covariants:
		cov_df = variant_df[(variant_df['Variant'] == cov)]
		isolates = cov_df[cov_df.columns[pd.Series(cov_df.columns).str.startswith('Isolates Count')]].iloc[:,:7].head(1)
		isolates = isolates.iloc[0].values.tolist()
		variants = cov_df[cov_df.columns[pd.Series(cov_df.columns).str.startswith('Variant Count')]].iloc[:,:7]
		variants = variants.fillna(0)
		variants.loc['Variant Total'] = variants.sum()
		count_sums = variants.loc['Variant Total'].values.tolist()
		cov_prevalence = [i / j for i, j in zip(count_sums, isolates)]
		prevalence_data = pd.DataFrame({'Variant': cov, prevalence_cols[5]: [cov_prevalence[5]], prevalence_cols[4]: [cov_prevalence[4]], prevalence_cols[3]: [cov_prevalence[3]], 
			prevalence_cols[2]: [cov_prevalence[2]], prevalence_cols[1]: [cov_prevalence[1]], prevalence_cols[0]: [cov_prevalence[0]]})
		growth_data = pd.DataFrame({'Variant': cov, growth_cols[5]: [div(cov_prevalence[5], cov_prevalence[6])], growth_cols[4]: [div(cov_prevalence[4],cov_prevalence[5])], 
			growth_cols[3]: [div(cov_prevalence[3],cov_prevalence[4])], growth_cols[2]: [div(cov_prevalence[2],cov_prevalence[3])], growth_cols[1]: [div(cov_prevalence[1],cov_prevalence[2])], 
			growth_cols[0]: [div(cov_prevalence[0],cov_prevalence[1])]})
		variant_data = pd.DataFrame({'Variant': cov, variant_cols[5]: [count_sums[5]], variant_cols[4]: [count_sums[4]], variant_cols[3]: [count_sums[3]], 
			variant_cols[2]: [count_sums[2]], variant_cols[1]: [count_sums[1]], variant_cols[0]: [count_sums[0]]})
		data_frames = [prevalence_data, growth_data, variant_data]
		data = reduce(lambda left, right: pd.merge(left, right, on = ['Variant'], how = 'outer'), data_frames)
		new_data.append(data)
	
	new_data = pd.concat(new_data, ignore_index = True)
	return(new_data)

# This function is meant to plot the covariants from either a PANGO lineage, WHO label, or inputted list.  
# Really, this function will plot the prevalence/growth rate trends of the covariants passed to the function.
# The covariants can either come from the results of one of the analysis options (substitution_ranking, functional_ranking,
# or composite ranking) or be user inputted.  For the sake of simplicity, the functional will only allow plotting maximum 10
# covariants at a time.  This function also offers the option to plot how the proportions of covariants that make up a lineage or a
# WHO clade change overtime.
def plot_covariants(variant_df, covariants, region, who = "", pango = "", name = {}, inputted = False):
	if (inputted):
		graph_type = input('Plot the change over time of covariant prevalence, growth rate, or proportion within the inputted list of covariants? [prevalence/growth]: ')
	else:
		graph_type = input('Plot the change over time of covariant prevalence, growth rate, or proportion within the PANGO/WHO clade? [prevalence/growth/proportion]: ')
	
	# Safety --- only allow potential plotting of covariants actually in the data file if supplying user list
	covariants = variant_df[(variant_df['Variant'].isin(covariants))]
	covariants = covariants.drop_duplicates(subset = ['Variant'])
	covariants = list(covariants['Variant'])

	if (pango):
		variant_df = variant_df[(variant_df['PANGO Lineage'] == pango)]
		variant_type = pango
	elif (who):
		variant_df = variant_df[(variant_df['WHO Label'] == who)]
		variant_df = plot_covariants_help(variant_df, covariants)
		variant_type = who
	elif (inputted == False):
		variant_df = variant_df[(variant_df['Variant'].isin(covariants))]
		variant_df = plot_covariants_help(variant_df, covariants)
		variant_type = ""
	else:
		variant_df = variant_df[(variant_df['Variant'].isin(covariants))]
		variant_df = plot_covariants_help(variant_df, covariants)
		variant_type = "User Inputted"


	if (graph_type == 'prevalence'):
		recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:6]
		recent_prevalence = recent_prevalence.fillna(0)
		covariant_prevalence_df = pd.concat([pd.DataFrame({'Variant': variant_df['Variant']}), recent_prevalence], axis = 1)
		covariant_prevalence_df = covariant_prevalence_df.reindex(sorted(covariant_prevalence_df.columns), axis = 1)
		covariant_prevalence_df = covariant_prevalence_df.set_index('Variant')
		covariant_prevalence_df['rowmax'] = covariant_prevalence_df.max(axis = 1)
		covariant_prevalence_df = covariant_prevalence_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		covariant_prevalence_df.columns = covariant_prevalence_df.columns.str.replace("Prevalence - ", "")
		x_val = np.asarray(covariant_prevalence_df.columns.tolist())
		indices = covariant_prevalence_df.index.tolist()
		if (name):
			label = [ name[k] for k in indices ]
		else:
			label = indices
		plt.stackplot(x_val, covariant_prevalence_df.loc[indices].values, labels = label)
		plt.suptitle('A Time Course of '+variant_type+' Covariate Prevalence Ratios from Sequences Isolated in '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Sequence Prevalence Ratio (Variant Count/Isolate Count)')
		plt.legend(loc = 'upper right', prop = {'size': 8}, bbox_to_anchor = (0.4, 1.05))
		plt.show()
	elif (graph_type == 'growth'):
		recent_growth_rate = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:6]
		recent_growth_rate = recent_growth_rate.fillna(0)
		covariant_growth_df = pd.concat([pd.DataFrame({'Variant': variant_df['Variant']}), recent_growth_rate], axis = 1)
		covariant_growth_df = covariant_growth_df.reindex(sorted(covariant_growth_df.columns), axis = 1)
		covariant_growth_df = covariant_growth_df.set_index('Variant')
		covariant_growth_df['rowmax'] = covariant_growth_df.max(axis = 1)
		covariant_growth_df = covariant_growth_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		covariant_growth_df.columns = covariant_growth_df.columns.str.replace("Growth - ", "")
		x_val = np.asarray(covariant_growth_df.columns.tolist())
		indices = covariant_growth_df.index.tolist()
		for i in indices:
			mask = np.isfinite(covariant_growth_df.loc[i].values)
			plt.plot(x_val[mask], covariant_growth_df.loc[i].values[mask], label =i)
		plt.suptitle('A Time Course of '+variant_type+' Covariate Growth Rates from Sequences Isolated in '+region, fontsize = 14)
		plt.xlabel('Time Intervals')
		plt.ylabel('Growth Rate')
		plt.legend(loc = 'upper right', prop = {'size': 6.5}, bbox_to_anchor = (0.4, 1.05))
		plt.show()
	elif (graph_type == 'proportion' and not inputted):
		recent_counts = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,:6]
		recent_counts = recent_counts.fillna(0)
		covariant_count_df = pd.concat([pd.DataFrame({'Variant': variant_df['Variant']}), recent_counts], axis = 1)
		covariant_count_df = covariant_count_df.reindex(sorted(covariant_count_df.columns), axis = 1)
		covariant_count_df = covariant_count_df.set_index('Variant')
		covariant_count_df.update(covariant_count_df.div(covariant_count_df.sum(axis=0),axis=1))
		covariant_count_df['rowmax'] = covariant_count_df.max(axis = 1)
		covariant_count_df = covariant_count_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		covariant_count_df.columns = covariant_count_df.columns.str.replace("Variant Count - ", "")
		x_val = np.asarray(covariant_count_df.columns.tolist())
		indices = covariant_count_df.index.tolist()
		plt.stackplot(x_val, covariant_count_df.loc[indices].values, labels = indices)
		plt.suptitle('A Time Course of Changes in Covariate Proportions that Make Up '+variant_type+' from Sequences Isolated In '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Covariant Proportional Abundance Relative to Clade')
		plt.legend(loc = 'upper right', prop = {'size': 6.5}, bbox_to_anchor = (0.4, 1.05))
		plt.show()
	else:
		sys.exit("Invalid graph type entry. Please try again and correctly specify the appropriate graph type (prevalence, growth, or proportion)")
		


# Plot the prevalence or growth rates of the user inputted lineage by country.  In this function the user
# will be able to specify the type of graph they want to display, which is either prevalence or growth rate change by
# country over time.  Note, this function could be supplied a pre-filtered variant_df dataframe that only includes
# data from certain countries, as potentially requested in main.py.  That said, the function will only work with data
# on countries in the supplied variant_df and report a graph for a maximum of 10 countires, the countries with the
# highest prevalence or growth rates for the variant.
def plot_single_lineage(variant_df, lineage):
	graph_type = input('Plot the changes overtime of the lineage prevalence or growth rate among the top/user inputted countries? [prevalence/growth]: ')
	variant_df = variant_df[variant_df[variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,0].name] > 10]
	variant_df = variant_df[(variant_df['PANGO Lineage'] == lineage)]

	if (graph_type == 'prevalence'):
		recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:6]
		recent_prevalence = recent_prevalence.fillna(0)
		lineage_prevalence_df = pd.concat([pd.DataFrame({'Country': variant_df['Country']}), recent_prevalence], axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.reindex(sorted(lineage_prevalence_df.columns), axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.set_index('Country')
		lineage_prevalence_df['rowmax'] = lineage_prevalence_df.max(axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		lineage_prevalence_df.columns = lineage_prevalence_df.columns.str.replace("Prevalence - ", "")
		x_val = np.asarray(lineage_prevalence_df.columns.tolist())
		indices = lineage_prevalence_df.index.tolist()
		for i in indices:
			mask = np.isfinite(lineage_prevalence_df.loc[i].values)
			plt.plot(x_val[mask], lineage_prevalence_df.loc[i].values[mask], label = i)
		plt.suptitle('A Time Course of the Proportion of Isolated Sequences with '+lineage+' by Country', fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Sequence Prevalence Ratio (Variant Count/Isolates Count)')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()
	elif (graph_type == 'growth'):
		recent_growth_rate = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:6]
		recent_growth_rate = recent_growth_rate.fillna(0)
		lineage_growth_df = pd.concat([pd.DataFrame({'Country': variant_df['Country']}), recent_growth_rate], axis = 1)
		lineage_growth_df = lineage_growth_df.reindex(sorted(lineage_growth_df.columns), axis = 1)
		lineage_growth_df = lineage_growth_df.set_index('Country')
		lineage_growth_df['rowmax'] = lineage_growth_df.max(axis = 1)
		lineage_growth_df = lineage_growth_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		lineage_growth_df.columns = lineage_growth_df.columns.str.replace("Growth Rate - ", "")
		x_val = np.asarray(lineage_growth_df.columns.tolist())
		indices = lineage_growth_df.index.tolist()
		for i in indices:
			mask = np.isfinite(lineage_growth_df.loc[i].values)
			plt.plot(x_val[mask], lineage_growth_df.loc[i].values[mask], label = i)
		plt.suptitle('Growth Rate Trends of '+lineage+' by Country', fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Growth Rate')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()
	else:
		sys.exit("Invalid graph type entry. Please try again and correctly specify the appropriate graph type (prevalence or growth)")


# Plot inputted PANGO lineage trends in prevalence or growth rates globally or for single country. In this function the user will 
# be able to specifiy the type of graph they want to display, which is either prevalence or growth rate.  The function 
# requires the 'World - Lineages' tab and a list of lineages to analyze. If
# this function is being called from a graph analysis option, then it will likely graph inputted lineages from the user
# for a global analysis or for a single country analysis.  Note, this function must take in a prefiltered variant_df to
# include only global prevalence/growth rates or prevalence/growth rates for a single country.  
def plot_lineages(variant_df, lineages, region):
	graph_type = input('Plot the changes over time of lineage prevalence or growth rate in '+region+'? [prevalence/growth]: ')
	variant_df = variant_df[(variant_df['PANGO Lineage'].isin(lineages))]

	if (graph_type == 'prevalence'):
		recent_prevalence = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Prevalence')]].iloc[:,:6]
		recent_prevalence = recent_prevalence.fillna(0)
		lineage_prevalence_df = pd.concat([pd.DataFrame({'Lineage': variant_df['PANGO Lineage']}), recent_prevalence], axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.reindex(sorted(lineage_prevalence_df.columns), axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.set_index('Lineage')
		lineage_prevalence_df['rowmax'] = lineage_prevalence_df.max(axis = 1)
		lineage_prevalence_df = lineage_prevalence_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(9)
		lineage_prevalence_df.loc['Other'] = 1 - lineage_prevalence_df.sum(axis = 0)
		lineage_prevalence_df.columns = lineage_prevalence_df.columns.str.replace("Prevalence - ", "")
		x_val = np.asarray(lineage_prevalence_df.columns.tolist())
		indices = lineage_prevalence_df.index.tolist()
		plt.stackplot(x_val, lineage_prevalence_df.loc[indices].values, labels = indices)
		plt.suptitle('A Time Course of Lineage Proportions from Sequence Isolates Extracted in '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Sequence Prevalence Ratio (Variant Count/Isolates Count)')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()

	elif (graph_type == 'growth'):
		recent_growth_rate = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,:6]
		recent_growth_rate = recent_growth_rate.fillna(0)
		lineage_growth_df = pd.concat([pd.DataFrame({'Lineage': variant_df['PANGO Lineage']}), recent_growth_rate], axis = 1)
		lineage_growth_df = lineage_growth_df.reindex(sorted(lineage_growth_df.columns), axis = 1)
		lineage_growth_df = lineage_growth_df.set_index('Lineage')
		lineage_growth_df['rowmax'] = lineage_growth_df.max(axis = 1)
		lineage_growth_df = lineage_growth_df.sort_values(by = 'rowmax', ascending = False).drop(labels = 'rowmax', axis = 1).head(10)
		lineage_growth_df.columns = lineage_growth_df.columns.str.replace("Growth Rate - ", "")
		x_val = np.asarray(lineage_growth_df.columns.tolist())
		indices = lineage_growth_df.index.tolist()
		for i in indices:
			mask = np.isfinite(lineage_growth_df.loc[i].values)
			plt.plot(x_val[mask], lineage_growth_df.loc[i].values[mask], label = i)
		plt.suptitle('Lineage Growth Rate Trends Over Time from Sequence Isolates Extracted in '+region, fontsize = 14)
		plt.xlabel('Year and Month')
		plt.ylabel('Growth Rate')
		plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.10), ncol = 5, fancybox = True, shadow = True)
		plt.show()
	#elif (graph_type == 'emergence'):
		#print("Emergence graph to come!")
		#for i in range(interval):
			#variant_df = variant_df[variant_df[variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Variant Count')]].iloc[:,i].name] > 10]
			#recent_growth_rates = variant_df[variant_df.columns[pd.Series(variant_df.columns).str.startswith('Growth Rate')]].iloc[:,i:i+interval]
			#headers = ["PANGO Lineage", "Growth Rates"]
			#data = []
			#for j in range(interval):
				#growth_data = pd.DataFrame({'PANGO Lineage': variant_df['PANGO Lineage'], 'Growth Rates': recent_growth_rates[recent_growth_rates.iloc[:,j].name]})
				#data.append(growth_data)
			#lineage_growth_data = pd.concat(data, ignore_index = True)
			#lineage_growth_data = lineage_growth_data[(lineage_growth_data['Growth Rates'] > 15)]
			#emergence_ranking = lineage_growth_data.groupby('PANGO Lineage').size().reset_index(name = 'Emergence Score').sort_values(by = 'Emergence Score', ascending = False)
			#emergence_ranking = emergence_ranking.reset_index(drop=True)

	else:
		sys.exit("Invalid graph type entry. Please try again and correctly specify the appropriate graph type (prevalence or growth)")
