# -*- coding: utf-8 -*-

import codecs
import ast
import argparse
import sys
import re
import pandas as pd
import numpy as np
import VariantAnalysis as va
import VariantPlots as vp


# Report the graph usage to the commandline
def graph_args_usage():
	print('\n')
	print("GRAPH USAGE: The following are permitted arguments for running the 'graph' analysis option")
	print('\n')
	print("(1) python main.py --filename [Emerging Variants Report] --analysis graph")
	print('---')
	print("(2) python main.py --filename [Emerging Variants Report] --analysis graph --lineage [PANGO Lineage]")
	print("---")
	print("(3) python main.py --filename [Emerging Variants Report] --analysis graph --lineage_file [Lineages TXT File]")
	print('---')
	print("(4) python main.py --filename [Emerging Variants Report] --analysis graph --country [Country]")
	print('---')
	print("(5) python main.py --filename [Emerging Variants Report] --analysis graph --lineage_file [Lineages TXT File] --country [Country]")
	print('---')
	print("(6) python main.py --filename [Emerging Variants Report] --analysis graph --lineage [PANGO Lineage] --country_file [Country TXT File]")
	print("---")
	print("(7) python main.py --filename [Emerging Variants Report] --analysis graph --[PANGO/WHO/covariants] [PANGO Lineage/WHO Label/Covariants TXT File]")
	print('---')
	print("(8) python main.py --filename [Emerging Variants Report] --analysis graph --[PANGO/WHO/covariants] [PANGO Lineage/WHO Label/Covariants TXT File] --country [Country]")
	print('---')
	print("(9) python main.py --filename [Emerging Variants Report] --analysis graph --mutations [Mutations TXT File]")
	print('---')
	print("(10) python main.py --filename [Emerging Variants Report] --analysis graph --domain [NTD/RBD/Spike/Other]")
	print('---')
	print("(11) python [Everything in (9) or (10)] --country [Country]")
	print('---')
	print("NOTE: Options 1-6 are for graphing PANGO Lineage trends, 7-9 is for graphing covariant trends, and 10-11 is for point mutation trends")
	print('\n')
	print("See ReadMe.md for details on graph analysis.")
	print('\n')

# Report the program usage to the command line
def program_usage():
	print('\n')
	print("RANKING USAGE: Commandline must have the minimum arguments: python main.py --filename [Emerging Variants Report] --analysis [Analysis Option]")
	print('\n')
	print("ANALYSIS OPTION: The following are arguments permitted after the '--analysis' argument")
	print('\n')
	print("(1) lineage_ranking: Rank PANGO lineages by their Emergence Score within entire Emerging Variants Report. These are 'growing' lineages")
	print("(a) python main.py --filename [Emerging Variants Report] --analysis lineage_ranking")
	print("(b) python [Everything in (a)] --interval [Interval]")
	print("(c) python [Everything in (a)] --country [Country]")
	print("(d) python [Everything in (a)] --interval [Interval] --country [Country]")
	print('---')
	print("(2) sequence_ranking: Rank covariants within a PANGO Lineage, WHO Label, inputted file, or entire Emerging Variants Report based on the Sequence Prevalence Score")
	print("(a) python main.py --filename [Emerging Variants Report] --analysis sequence_ranking")
	print("(b) python main.py --filename [Emerging Variants Report] --analysis sequence_ranking --[PANGO/WHO/covariants] [PANGO Lineage/WHO Label/Covariants TXT File]")
	print("(c) python [Everything in (a) or (b)] --interval [Interval]")
	print("(d) python [Everything in (a) or (b)] --country [Country]")
	print("(e) python [Everything in (a) or (b)] --interval [Interval] --country [Country]")
	print('---')
	print("(3) functional_ranking: Rank covariants within a PANGO Lineage, WHO Label, inputed file, or entire Emerging Variants Report based on the Predicted Functional Impact Score")
	print("(a) python main.py --filename [Emerging Variants Report] --analysis functional_ranking")
	print("(b) python main.py --filename [Emerging Variants Report] --analysis functional_ranking --[PANGO/WHO/covariants] [PANGO Lineage/WHO Label/Covariants TXT File]")
	print("(c) python [Everything in (a) or (b)] --interval [Interval]")
	print("(d) python [Everything in (a) or (b)] --country [Country]")
	print("(e) python [Everything in (a) or (b)] --interval [Interval] --country [Country]")
	print('---')
	print("(4) composite_ranking: Rank covariants within a PANGO Lineage, WHO Label, inputted file, or Entire Emerging Variants Report based on the Composite Score")
	print("(a) python main.py --filename [Emerging Variants Report] --analysis composite_ranking")
	print("(b) python main.py --filename [Emerging Variants Report] --analysis composite_ranking --[PANGO/WHO/covariants] [PANGO Lineage/WHO Label/Covariants TXT File]")
	print("(c) python [Everything in (a)] --interval [Interval]")
	print("(d) python [Everything in (a)] --country [Country]")
	print("(e) python [Everything in (a)] --interval [Interval] --country [Country]")
	print('---')
	print("(5) mutation_ranking: Rank AA substitutions within the NTD, RBD, outside of those two domains, or witin the entire Spike based on the Sequence Prevalence Score")
	print("(a) python main.py --filename [Emerging Variants Report] --analysis mutation_ranking --domain [NTD/RBD/Spike/Other]")
	print("(b) python main.py --filename [Emerging Variants Report] --analysis mutation_ranking --mutations [Mutations TXT File]")
	print("(b) python [Everything in (a) or (b)] --interval [Interval]")
	print("(c) python [Everything in (a) or (b)] --country [Country]")
	print("(e) python [Everything in (a) or (b)] --[interval/domain/mutations/country] [Interval/Country/[Domains]/Mutations TXT File")
	print('---')
	print("(6) graph: Graph trends of PANGO Lineages, trends of covariants within a PANGO Lineage, WHO Label, or inputed file, or trends of single point mutations")
	print('---')
	print("(7) help: Print out the GENERAL USAGE and GRAPH USAGE")
	print('\n')
	print("See ReadMe.md for more specific details regarding usage.")
	print('\n')


# Used for opening the covariants file.  This file must have at least one column named 'Variant'
# WARNING: Must makes sure entries in Variant Column match the variants in the spreadsheet. Otherwise, analysis results will be off.
def open_covariants_file(file): 
	try:
		covariants_file = pd.read_csv(file, sep = '\t', encoding = 'latin-1')
		if ('Variant' in covariants_file.columns):
			if ('Name' in covariants_file.columns):
				covariants_file['Name'] = covariants_file['Name'].str.encode('ascii', 'ignore').str.decode('ascii')
				covariants_file['Name'] = covariants_file['Name'].str.replace('Ê', '')
			covariants_file['Variant'] = covariants_file['Variant'].str.replace(' ', '')
			covariants_file['Variant'] = covariants_file['Variant'].str.replace('Ê', '')
			return(covariants_file)
		else:
			sys.exit("Invalid covariants file format. Inputted covariants file must be tabular with one column labeled 'Variant'.  Optional 'Name' or 'Identifier' columns permitted. See ReadMe.md for more details.")
	except:
		raise Exception("Could not open covariants file")

# Used for opening the file with list of countries, lineages, and single point mutations
def open_optional_file(file):
	try:
		entries = open(file).read().splitlines()
		return(entries)
	except:
		raise Exception("Could not open optional TXT file")

# Opens the Emerging Variants Report and returns the data from the spreadsheet
def open_emerging_variants(file):
	try:
		print("Opening the Emerging Variants Report ...")
		variants_xlsx = pd.ExcelFile(file)
		return(variants_xlsx)
	except:
		raise Exception("Could not open file for Emerging Variants Report")



if __name__ == '__main__':
	analysis_options = ['lineage_ranking', 'sequence_ranking', 'functional_ranking', 'composite_ranking', 'mutation_ranking', 'graph', 'help']
	domains = ['NTD', 'RBD', 'Spike', 'Other']
	
	if (len(sys.argv) < 4 or len(sys.argv) > 10):
		sys.exit(program_usage())
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--filename', dest = 'filename', type = str)
	parser.add_argument('--analysis', dest = 'analysis', type = str)
	parser.add_argument('--PANGO', dest = 'pango', type = str)
	parser.add_argument('--WHO', dest = 'who', type = str)
	parser.add_argument('--covariants', dest = 'covariant', type = str)
	parser.add_argument('--interval', dest = 'interval', type = int)
	parser.add_argument('--lineage_file', dest = 'lineage_file', type = str)
	parser.add_argument('--lineage', dest = 'lineage', type = str)
	parser.add_argument('--country_file', dest = 'country_file', type = str)
	parser.add_argument('--country', dest = 'country', type = str)
	parser.add_argument('--domain', dest = 'domain', type = str)
	parser.add_argument('--mutations', dest = 'mutation', type = str)
	args = parser.parse_args()
	
	# Ensure proper arguments to the commandline
	if (not args.filename):
		sys.exit(program_usage())
	if (args.analysis not in analysis_options):
		sys.exit(program_usage())
	if (args.analysis == 'lineage_ranking' and (args.pango or args.who or args.covariant)):
		sys.exit(program_usage())
	if (args.pango and args.who and args.covariant):
		sys.exit(program_usage())
	if (args.pango and (args.who or args.covariant)):
		sys.exit(program_usage())
	if (args.who and (args.pango or args.covariant)):
		sys.exit(program_usage())
	if (args.covariant and (args.pango or args.who)):
		sys.exit(program_usage())
	if ((args.analysis != 'mutation_ranking' and args.analysis != 'graph') and (args.domain or args.mutation)):
		sys.exit(program_usage())
	if (args.analysis == 'mutation_ranking' and (not args.domain and not args.mutation)):
		sys.exit(program_usage())
	if (args.analysis == 'mutation_ranking' and (not args.mutation) and (args.domain not in domains)):
		sys.exit(program_usage())
	if (args.analysis == 'help'):
		program_usage()
		graph_args_usage()
		sys.exit()
	if (args.analysis != 'graph' and (args.lineage_file or args.lineage or args.country_file)):
		sys.exit(program_usage())
	if (args.analysis == 'graph' and (args.lineage_file and args.country_file)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and (args.lineage and args.country)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and (not args.lineage and args.country_file)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and args.lineage and (args.pango or args.who or args.covariant)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and (args.lineage or args.lineage_file or args.country_file) and (args.domain or args.mutation)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and (args.pango or args.who or args.covariant) and (args.domain or args.mutation)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and (args.domain and args.mutation)):
		sys.exit(graph_args_usage())
	if (args.analysis == 'graph' and args.domain and (args.domain not in domains)):
		sys.exit(graph_args_usage())
	if (args.interval):
		interval = args.interval
	else:
		interval = 4

	
	# Filter the data by country or countries if needed
	variants_xlsx = open_emerging_variants(args.filename)
	if (not args.country):
		variants = pd.read_excel(variants_xlsx, 'World - Variants')
		lineages = pd.read_excel(variants_xlsx, 'World - Lineages')
		mutations = pd.read_excel(variants_xlsx, 'AA Mutations')
		analysis_variants = variants[(variants['Country'] != 'All') & (variants['Country'] != 'Unknown')]
		graph_variants = variants[(variants['Country'] == 'All')]
		lineages = lineages[(lineages['Country'] == 'All')]
		analysis_subs = mutations[(mutations['Country'] != 'All') & (mutations['Country'] != 'Unknown')]
		graph_subs = mutations[(mutations['Country'] == 'All')]
		region = "World"
	elif (args.country == "USA"):
		analysis_variants = pd.read_excel(variants_xlsx, 'USA - Variants')
		lineages = pd.read_excel(variants_xlsx, 'USA - Lineages')
		mutations = pd.read_excel(variants_xlsx, 'AA Mutations')
		graph_variants = analysis_variants[(analysis_variants['Region'] == 'All')]
		lineages = lineages[(lineages['Region'] == 'All')]
		graph_subs = mutations[(mutations['Country'] == 'USA')]
		region = "USA"
	else:
		variants = pd.read_excel(variants_xlsx, 'World - Variants')
		lineages = pd.read_excel(variants_xlsx, 'World - Lineages')
		mutations = pd.read_excel(variants_xlsx, 'AA Mutations')
		analysis_variants = variants[(variants['Country'] == args.country)]
		graph_variants = analysis_variants
		lineages = lineages[(lineages['Country'] == args.country)]
		analysis_subs = mutations[(mutations['Country'] == args.country)]
		graph_subs = analysis_subs
		region = args.country

	if (args.country_file):
		countries = open_optional_file(args.country_file)
		lineages = pd.read_excel(variants_xlsx, 'World - Lineages')
		lineages = lineages[(lineages['Country'].isin(countries))]


	# The following code conducts the range of analysis options
	if (args.analysis == 'lineage_ranking'):
		emergence_ranking = va.lineage_ranking(analysis_variants, interval)
		emergence_ranking.to_csv('lineage_ranking_'+region+'.tsv', sep = '\t', index = False)
		print(emergence_ranking)
		print('\n')
		graph = input('Display a graph of trends over time for top emerging lineages in '+region+'? [y/n]: ')
		if (graph == 'y'):
			vp.plot_lineages(lineages, emergence_ranking['PANGO Lineage'], region)
	
	elif (args.analysis == 'sequence_ranking'):
		if (args.covariant):
			covariants_file = open_covariants_file(args.covariant)
			covariants = list(covariants_file['Variant'])
			covariant_score = va.sequence_ranking(analysis_variants, interval, covariants = covariants)
			covariant_score = pd.merge(covariants_file, covariant_score, on = "Variant")
			covariant_score = covariant_score.sort_values(by = 'Sequence Prevalence Score', ascending = False).reset_index(drop = True)
			covariant_score.to_csv('inputted_covariants_sequence_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariant_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top user inputted covariants from the sequence ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				if ('Name' in covariants_file.columns):
					names = list(covariants_file['Name'])
					name_dict = dict(zip(covariants, names))
					vp.plot_covariants(graph_variants, covariants, region, name = name_dict)
				else:
					vp.plot_covariants(graph_variants, covariants, region)
		elif (args.pango):
			pango_lineage = args.pango
			if (analysis_variants['PANGO Lineage'].isin([pango_lineage]).any()):
				covariant_score = va.sequence_ranking(analysis_variants, interval, lineage = [pango_lineage])
				covariant_score.to_csv(pango_lineage+'_sequence_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the sequence ranking for '+pango_lineage+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, pango = pango_lineage)
			else:
				sys.exit("PANGO Lineage is invalid and not found. Please try again an input a valid PANGO Lineage.")
		elif (args.who):
			who_label = args.who
			if (analysis_variants['WHO Label'].isin([who_label]).any()):
				covariant_score = va.sequence_ranking(analysis_variants, interval, who = [who_label])
				covariant_score.to_csv(who_label+'_sequence_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the substitution ranking for '+who_label+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, who = who_label)
			else:
				sys.exit("WHO Label is invalid and not found. Please try again and input a valid WHO Label.")
		else:
			covariants_score = va.sequence_ranking(analysis_variants, interval)
			covariants_score.to_csv('emerging_covariants_sequence_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariants_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top covariants from the sequence ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				vp.plot_covariants(graph_variants, covariants_score['Variant'].tolist(), region, inputted = True)
	
	elif (args.analysis == 'functional_ranking'):
		sfocs = pd.read_excel(variants_xlsx, 'SFoCs')
		sfocs = sfocs[(sfocs['Protein'] == 'Spike')]
		if (args.covariant):
			covariants_file = open_covariants_file(args.covariant)
			covariants = list(covariants_file['Variant'])
			covariant_score = va.functional_ranking(analysis_variants, sfocs, covariants = covariants)
			covariant_score = pd.merge(covariants_file, covariant_score, on = "Variant")
			covariant_score = covariant_score.drop_duplicates(subset = ['Variant']).sort_values(by = 'Functional Impact Score', ascending = False).reset_index(drop = True)
			covariant_score.to_csv('inputted_covariants_functional_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariant_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top inputted covariants from the functional_ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				if ('Name' in covariants_file.columns):
					names = list(covariants_file['Name'])
					name_dict = dict(zip(covariants, names))
					vp.plot_covariants(graph_variants, covariants, region, name = name_dict)
				else:
					vp.plot_covariants(graph_variants, covariants, region)
		elif (args.pango):
			pango_lineage = args.pango
			if (analysis_variants['PANGO Lineage'].isin([pango_lineage]).any()):
				covariant_score = va.functional_ranking(analysis_variants, sfocs, lineage = [pango_lineage])
				covariant_score.to_csv(pango_lineage+'_functional_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the functional_ranking for '+pango_lineage+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, pango = pango_lineage)
			else:
				sys.exit("PANGO Lineage is invalid and not found. Please try again an input a valid PANGO Lineage.")
		elif (args.who):
			who_label = args.who
			if (analysis_variants['WHO Label'].isin([who_label]).any()):
				covariant_score = va.functional_ranking(analysis_variants, sfocs, who = [who_label])
				covariant_score.to_csv(who_label+'_functional_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the functional_ranking for '+who_label+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, who = who_label)
			else:
				sys.exit("WHO Label is invalid and not found. Please try again and input a valid WHO Label.")
		else:
			print('WARNING! This computation may take some time! Reconsider if necessary. Control + C to abort.')
			analysis_variants = analysis_variants[analysis_variants[analysis_variants[analysis_variants.columns[pd.Series(analysis_variants.columns).str.startswith('Variant Count')]].iloc[:,0].name] > 10]
			covariants_score = va.functional_ranking(analysis_variants, sfocs)
			covariants_score.to_csv('emerging_covariants_functional_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariants_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top covariants from the functional ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				vp.plot_covariants(graph_variants, covariants_score['Variant'].tolist(), region, inputted = True)
	
	elif (args.analysis == 'composite_ranking'):
		sfocs = pd.read_excel(variants_xlsx, 'SFoCs')
		sfocs = sfocs[(sfocs['Protein'] == 'Spike')]
		if (args.covariant):
			covariants_file = open_covariants_file(args.covariant)
			covariants = list(covariants_file['Variant'])
			covariant_score = va.composite_ranking(analysis_variants, sfocs, interval, covariants = covariants)
			covariant_score = pd.merge(covariants_file, covariant_score, on = "Variant")
			covariant_score = covariant_score.drop_duplicates(subset = ['Variant']).sort_values(by = 'Composite Score', ascending = False).reset_index(drop = True)
			covariant_score.to_csv('inputted_covariants_composite_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariant_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top inputted covariants from the composite_ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				if ('Name' in covariants_file.columns):
					names = list(covariants_file['Name'])
					name_dict = dict(zip(covariants, names))
					vp.plot_covariants(graph_variants, covariants, region, name = name_dict, inputted = True)
				else:
					vp.plot_covariants(graph_variants, covariants, region, inputted = True)
		elif (args.pango):
			pango_lineage = args.pango
			if (analysis_variants['PANGO Lineage'].isin([pango_lineage]).any()):
				covariant_score = va.composite_ranking(analysis_variants, sfocs, interval, lineage = [pango_lineage])
				covariant_score.to_csv(pango_lineage+'_composite_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the composite_ranking for '+pango_lineage+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, pango = pango_lineage)
			else:
				sys.exit("PANGO Lineage is invalid and not found. Please try again an input a valid PANGO Lineage.")
		elif (args.who):
			who_label = args.who
			if (analysis_variants['WHO Label'].isin([who_label]).any()):
				covariant_score = va.composite_ranking(analysis_variants, sfocs, interval, who = [who_label])
				covariant_score.to_csv(who_label+'_composite_ranking_'+region+'.tsv', sep = '\t', index = False)
				print(covariant_score)
				print('\n')
				graph = input('Display a graph of trends over time for the top covariants from the composite_ranking for '+who_label+' in '+region+'? [y/n]: ')
				if (graph == 'y'):
					vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, who = who_label)
			else:
				sys.exit("WHO Label is invalid and not found. Please try again and input a valid WHO Label.")
		else:
			covariants_score = va.composite_ranking(analysis_variants, sfocs, interval)
			covariants_score.to_csv('emerging_covariants_composite_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(covariants_score)
			print('\n')
			graph = input('Display a graph of trends over time for the top covariants from the composite ranking in '+region+'? [y/n]: ')
			if (graph == 'y'):
				vp.plot_covariants(graph_variants, covariants_score['Variant'].tolist(), region, inputted = True)

	elif (args.analysis == 'mutation_ranking'):
		domain = args.domain
		if (args.mutation):
			mutation_file = open_optional_file(args.mutation)
			mutation_ranking = va.mutation_ranking(analysis_subs, interval, domain, mutation_file)
			mutation_ranking.to_csv('inputted_mutations_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(mutation_ranking)
		else:
			mutation_ranking = va.mutation_ranking(analysis_subs, interval, domain)
			mutation_ranking.to_csv(domain+'_mutations_ranking_'+region+'.tsv', sep = '\t', index = False)
			print(mutation_ranking)
		graph = input('Display a graph of trends over time for the top point mutations ranked within '+domain+' in '+region+'? [y/n]: ')
		if (graph == 'y'):
			vp.plot_mutations(graph_subs, region, domain, mutation_ranking['Variant'].tolist())


	elif (args.analysis == 'graph'):
		if (args.lineage_file):
			pango_lineages = open_optional_file(args.lineage_file)
			if ((not args.country) and (len(pango_lineages) == 1)):
				lineages = pd.read_excel(variants_xlsx, 'World - Lineages')
				lineages = lineages[(lineages['Country'] != 'All') & (lineages['Country'] != 'Unknown')]
				vp.plot_single_lineage(lineages, pango_lineages[0])
			else:
				vp.plot_lineages(lineages, pango_lineages, region)
		elif (args.lineage):
			pango_lineage = args.lineage
			if (not args.country_file):
				lineages = pd.read_excel(variants_xlsx, 'World - Lineages')
				lineages = lineages[(lineages['Country'] != 'All') & (lineages['Country'] != 'Unknown')]
				vp.plot_single_lineage(lineages, pango_lineage)
			else:
				vp.plot_single_lineage(lineages, pango_lineage)
		elif (args.covariant):
			covariants_file = open_covariants_file(args.covariant)
			covariants = list(covariants_file['Variant'])
			if (analysis_variants['Variant'].isin(covariants).any()):
				if ('Name' in covariants_file.columns):
					names = list(covariants_file['Name'])
					name_dict = dict(zip(covariants, names))
					vp.plot_covariants(graph_variants, covariants, region, name = name_dict, inputted = True)
				else:
					vp.plot_covariants(graph_variants, covariants, region, inputted = True)
			else:
				sys.exit("All covariants are invalid and not found. Please try again and input valid covariants.")
		elif (args.pango):
			pango_lineage = args.pango
			if (analysis_variants['PANGO Lineage'].isin([pango_lineage]).any()):
				covariant_score = va.sequence_ranking(analysis_variants, interval, lineage = [pango_lineage])
				vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, pango = pango_lineage, inputted = False)
			else:
				sys.exit("PANGO Lineage is invalid and not found. Please try again and input a valid PANGO Lineage.")
		elif (args.who):
			who_label = args.who
			if (analysis_variants['WHO Label'].isin([who_label]).any()):
				covariant_score = va.sequence_ranking(analysis_variants, interval, who = [who_label])
				vp.plot_covariants(graph_variants, covariant_score['Variant'].tolist(), region, who = who_label, inputted = False)
			else:
				sys.exit("WHO Label is invalid and not found. Please try again and input a valid WHO Label.")
		elif (args.mutation):
			mutations = open_optional_file(args.mutation)
			if (graph_subs['Variant'].isin(mutations).any()):
				vp.plot_mutations(graph_subs, region, domain = "User Inputted", mutations = mutations)
			else:
				sys.exit("Point mutations are invalid and not found. Please try again and input a valid list of SARS-CoV-2 Spike mutations.")
		elif (args.domain):
			domain = args.domain
			vp.plot_mutations(graph_subs, region, domain)
		else:
			recent_prevalence = analysis_variants[analysis_variants.columns[pd.Series(analysis_variants.columns).str.startswith('Prevalence')]].iloc[:,:6]
			recent_prevalence = recent_prevalence.fillna(0)
			lineage_prevalence_df = pd.concat([pd.DataFrame({'Lineage': analysis_variants['PANGO Lineage']}), recent_prevalence], axis = 1)
			lineage_prevalence_df = lineage_prevalence_df.reindex(sorted(lineage_prevalence_df.columns), axis = 1).drop_duplicates(subset = ['Lineage'])
			pango_lineages = lineage_prevalence_df['Lineage']
			vp.plot_lineages(lineages, pango_lineages, region)




