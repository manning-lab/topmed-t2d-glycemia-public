import argparse
parser = argparse.ArgumentParser(description='Get a list of sample IDs for all input cohorts within a firecloud data model and print to a text file, one sample id per line',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--cohorts', nargs=1, default='', help='common separated list of cohorts that you would like the sample ids from')
parser.add_argument('--outfile_pref', nargs=1, default='', help='label for output file of sample ids')

args = parser.parse_args()
args.cohorts = args.cohorts[0].split(',')
out_file = args.outfile_pref[0]+'.txt'


from firecloud import fiss
samples = fiss.fapi.get_entities('topmed-shared','topmed-shared', 'sample').json()
sample_study_gen = (s['attributes']['participant']['entityName'] for s in samples if s['attributes']['study'] in cohorts)

with open(out_file, 'w') as f:
	for p in sample_study_gen:
		f.write(p+'\n')