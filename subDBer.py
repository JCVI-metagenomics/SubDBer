#!/usr/local/bin/python


## Dependencies ##

import sys
import re
import os
import tarfile
import urllib
import subprocess


## Help information ##

print '\n'.join((
	'',
	'subDBer.py',
	'Author: Qiyun Zhu (qzhu@jcvi.org)',
	'Affiliation: J. Craig Venter Institute',
	'Last updated: Oct. 4, 2015',
	''
	))

if len(sys.argv) == 1:
	print '\n'.join((
		'Usage:',
		'  python subDBer.py [-in dbname] [-out dbname] [-outfmt fasta/blast/gilist/alias] [-title dbtitle] [-within taxids] [-exclude taxids] [-rank rank] [-size number] [-keep taxids] [-prefer auto] [-quiet] [-help]',
		'',
		'Example:',
		'  python subDBer.py -in nt -out sub_nt -outfmt blast -within 2 -exclude 1117 -rank genus -size 1 -keep 816,838,1263',
		'',
		'Explanation:',
		'  This command takes the NCBI nt database as input -> starts with all organisms from Bacteria (TaxID: 2) except for Cyanobacteria (1117) -> picks one representative organism per genus -> except for Bacteroides (816), Prevotella (838) and Ruminococcus (1263), in which all organisms are included -> creates a new BLAST database sub_nt that contains sequence data from the selected organisms.',
		'',
		'Use "-help" to display details of program parameters.',
		''
		))
	sys.exit()

helpinfo = '\n'.join((
	'Details of parameters:',
	'  -in <input database name>',
	'    The input database should be a compiled BLAST database. You may download a database of choice from the NCBI FTP site, or create one based on custom sequences.',
	'    The database may be defined by its name (e.g., env_nr), if it is already stored in the local file system and included in the environmental variable. Otherwise, you may enter the path to the database (e.g., /home/me/db/env_nr).',
	'    If the database is not found in the local file system but is available from the NCBI server (e.g., nt), SubDBer will download the latest copy.',
	'  -out <output database name>',
	'    Name of the output database. If omitted, SubDBer will add "_sub" after the name of the input database',
	'  -outfmt <format of output database>',
	'    Options:',
	'    1) fasta (default): a multisequence sequence FASTA file',
	'    2) blast: a compiled BLAST database',
	'    3) alias: an alias file (.nal) that represents a subset of the original BLAST database. It can be used to replace the "db" parameter in future BLAST searches. This saves the need for extra disk space to store the actual subsampled BLAST database.',
	'    4) gilist: a list of GI\'s (.gil). It can be used as the "gilist" parameter of future BLAST searches. Similar to alias.',
	'  -title <title of output database>',
	'    If omitted, the title will be the same as the input database (if any), or "untitled"',
	'  -within <TaxIDs to start with, separated by comma>',
	'    All child taxa under the designated TaxIDs will be considered as the initial population for subsampling.',
	'    TaxIDs can be looked up from the NCBI taxonomy database.',
	'  -exclude <TaxIDs to exlude>',
	'    All child taxa will be excluded from the initial population.',
	'  -keep <TaxIDs to keep>',
	'    All child taxa will be retained (not subject to subsampling).',
	'  -rank <taxonomic rank as unit for subsampling>',
	'    Options: phylum, class, order, family, genus, species and other ranks that are acceptable in the NCBI taxonomy database.',
	'    If omitted, all taxa will be retained (no subsampling).',
	'  -size <number of taxa>',
	'    Up to this number of taxa will be sampled from each taxonomic group as defined by -rank. Default: 1.',
	'    Note: this sampling process is not random. The rule is that taxa with more synonyms recorded in the database have higher priority to be sampled.',
	'  -prefer <preferred TaxIDs>',
	'    A file containing one TaxID per line. These organisms will be preferred during the subsampling process.',
	'    If the value is "auto", the program will refer to the standard representative genome list at the NCBI server.',
	'  -quiet',
	'    SubDBer will complete all procedures automatically (non-interactively).',
	'  -help',
	'    Print a detailed list of all parameters.',
	''
))


## Parameters ##

(indb, outdb) = ('', '')
(rank, size) = ('', 1)
(within, exclude, keep, prefer) = ([], [], [], '')
(title, dbtype) = ('untitled', 'nucl')
outfmt = 'fasta'
quiet = 0


## Global variables ##

taxdumps = {}
	# the NCBI taxonomy database
	# id: name, parent, rank, children, names
taxa = {}
	# all taxa mentioned in the input database
	# id: seqs, length, selected
recs = {}
	# all records stored in the input database
	# oid: gi, accn, taxid, length
groups = {}
	# taxonomic groups on the designated rank in the input database
	# id: number of taxa


## Functions ##

def iterate_filter(id, val):
	if taxdumps.has_key(id):
		taxdumps[id]['in'] = val
		if taxdumps[id]['children']:
			for cid in taxdumps[id]['children']:
				iterate_filter(cid, val)


## Read parameters ##

args = len(sys.argv)
for i in range(1, args):
	if sys.argv[i][0] != '-': continue
	if args > i+1 and sys.argv[i+1][0] != '-':
		(arg, val) = (sys.argv[i], sys.argv[i+1])
		if arg == '-in': indb = val
		elif arg == '-out': outdb = val
		elif arg == '-outfmt': outfmt = val
		elif arg == '-title': title = val
		elif arg == '-within': within = val.split(',')
		elif arg == '-exclude': exclude = val.split(',')
		elif arg == '-keep': keep = val.split(',')
		elif arg == '-prefer': prefer = val
		elif arg == '-rank': rank = val
		elif arg == '-size': size = int(val)
	else:
		arg = sys.argv[i]
		if arg == '-quiet': quiet = 1
		elif arg == '-help': sys.exit(helpinfo)

if indb:
	if not outdb: outdb = indb[indb.rfind('/')+1:] + '_sub'
else:
	sys.exit('Please specify an input BLAST database.')


## Step 1: Download the latest NCBI taxonomy database ##

if os.path.isfile('names.dmp') and os.path.isfile('nodes.dmp'):
	print 'The taxonomy database is present in the current directory. subDBer will use it.'
	if not quiet: raw_input('Press any key to continue...')
else:
	print 'The taxonomy database does not exist or is incomplete. subDBer will download it from the NCBI server.'
	if not quiet: raw_input('Press any key to continue...')
	for file in ['names.dmp', 'nodes.dmp', 'taxdump.tar.gz']:
		if os.path.isfile(file): os.remove(file)
	print '  Downloading...'
	urllib.urlretrieve ('ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz', 'taxdump.tar.gz')
	print '  Extracting...'
	tar = tarfile.open('taxdump.tar.gz', 'r:gz')
	tar.extract('names.dmp', '.')
	tar.extract('nodes.dmp', '.')
	tar.close()
	os.remove('taxdump.tar.gz')
	if os.path.isfile('names.dmp') and os.path.isfile('nodes.dmp'):
		print 'The database is successfully downloaded and extracted.'
		if not quiet: raw_input('Press any key to continue...')
	else:
		sys.exit('Error: Failed to retrieve the database.')


## Step 2: Read the NCBI taxonomy database ##

print 'Reading taxonomy database...'
fnodes = open('nodes.dmp', 'r')
for line in fnodes:
	l = re.split('\s+\|\s+', line)
	taxdumps[l[0]] = {'parent': l[1], 'rank': l[2], 'children': [], 'names': 0}
fnodes.close()
for id in taxdumps:
	if (taxdumps[id]['parent'] and taxdumps.has_key(taxdumps[id]['parent'])):
		taxdumps[taxdumps[id]['parent']]['children'].append(id)
fnames = open('names.dmp', 'r')
for line in fnames:
	l = re.split('\s+\|\s+', line)
	if not taxdumps.has_key(l[0]): continue
	taxdumps[l[0]]['names'] += 1
	if re.search('scientific name\s*\|', line): taxdumps[l[0]]['name'] = l[1]
fnames.close()
print '  ' + str(len(taxdumps)) + ' taxa read.'

print 'Filtering taxonomy database...'
if within:
	within_names = []
	for id in within:
		if taxdumps.has_key(id):
			within_names.append(taxdumps[id]['name'] + ' (' + id + ')')
			iterate_filter(id, 1)
	for id in taxdumps.keys():
		if not taxdumps[id].has_key('in'):
			del taxdumps[id]
	print '  ' + str(len(taxdumps)) + ' taxa under ' + ', '.join(within_names) + ' are retained.'

if exclude:
	exclude_names = []
	excluded = 0
	for id in exclude:
		if taxdumps.has_key(id):
			exclude_names.append(taxdumps[id]['name'] + ' (' + id + ')')
			iterate_filter(id, 0)
	for id in taxdumps.keys():
		if taxdumps[id].has_key('in') and taxdumps[id]['in'] == 0:
			del taxdumps[id]
			excluded += 1
	print '  ' + str(excluded) + ' taxa under ' + ', '.join(exclude_names) + ' are excluded, remaining ' + str(len(taxdumps)) + ' taxa.'


## Step 3: Download the representative genome list ##

lprefer = {}
if prefer == 'auto':
	print 'SubDBer will refer to the NCBI representative genome list for preferred genomes.'
	preferfile = 'representative_genomes.txt'
	if os.path.isfile(preferfile):
		print 'The representative genomes list is present in the current directory. subDBer will use it.'
		if not quiet: raw_input('Press any key to continue...')
	else:
		print 'Downloading representative genome list...',
		sys.stdout.flush()
		try:
			urllib.urlretrieve ('http://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=refgenomes&download=on', preferfile)
			print 'done.'
		except:
			print 'failed.'
	if os.path.isfile(preferfile):
		print 'Reading representative genome list...'
		name2id = {}
		for i in taxdumps:
			if not taxdumps[i].has_key('name'): continue
			name = taxdumps[i]['name']
			if not name2id.has_key(name):
				name2id[name] = i
		fprefer = open(preferfile, 'r')
		for line in fprefer:
			if line[0] == '#': continue
			line = line.rstrip('\r\n')
			if not line: continue
			name = line.split('\t')[2]
			if name2id.has_key(name):
				lprefer[name2id[name]] = 1
		fprefer.close()
		print '  ' + str(len(lprefer)) + ' representative genomes read.'
elif prefer:
	if not os.path.isfile(prefer):
		sys.exit('Error: The preferred genome list does not exist.')
	print 'Reading preferred genome list...'
	fprefer = open(prefer, 'r')
	for line in fprefer:
		line = line.rstrip('\r\n')
		lprefer[line] = 1
	fprefer.close()
	print '  ' + str(len(lprefer)) + ' preferred genomes read.'


## Step 4: Read a designated BLAST database ##

print "Verifying input BLAST database..."
p = subprocess.Popen('blastdbcmd -db ' + indb + ' -info', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
out = p.stdout.read()
if 'command not found' in out: sys.exit('Error: blastdbcmd not available. Check if you have ncbi-blast+ correctly installed.')
if 'BLAST Database error: No alias or index' in out:
	print 'Error: The input BLAST database ' + indb + ' is missing or invalid.'
	p = subprocess.Popen('update_blastdb.pl --showall', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	if 'command not found' in out: sys.exit('Error: update_blastdb.pl not available. Check if you have ncbi-blast+ correctly installed.')
	ncbi_dbs = []
	for line in out.split('\n'):
		if line[:9] == 'Connected': continue
		ncbi_dbs.append(line.rstrip())
	indb = indb[indb.rfind('/')+1:]
	if not indb in ncbi_dbs: sys.exit('Error: Database ' + indb + ' is not available at the NCBI server either.')
	print 'Database ' + indb + ' is available from the NCBI server. SubDBer will download it to the current directory.'
	if not quiet: raw_input('Press any key to continue...')
	p = subprocess.call('update_blastdb.pl --decompress ' + indb, shell=True)
	print "Verifying downloaded BLAST database..."
	p = subprocess.Popen('blastdbcmd -db ' + indb + ' -info', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	if 'BLAST Database error: No alias or index' in out:
		print 'Error: The downloaded database ' + indb + ' is still missing or invalid.'
if 'total bases' in out: dbtype= 'nucl'
elif 'total residues' in out: dbtype= 'prot'
for line in out.split('\n'):
	if not line: continue
	if line[0] == '\t': line = '  ' + line[1:]
	print '  ' + line
	if line[:10] == 'Database: ' and title == 'untitled': title = line[10:]
p.stdout.close()
if ' ' in title: title = '"' + title + '"'

print "Reading input BLAST database..."
allrecs = 0
p = subprocess.Popen('blastdbcmd -db ' + indb + ' -entry all -outfmt \"%o %g %a %T %l\"', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
for line in p.stdout:
	line = line.rstrip('\r\n')
	if not line: continue
	allrecs += 1
	l = line.split()
	if not l[3]: continue
	if not taxdumps.has_key(l[3]): continue
		# only taxa recorded in the taxonomy database are retained.
	recs[l[0]] = {'gi': l[1], 'accn': l[2], 'taxid': l[3], 'length': l[4]}
	if taxa.has_key(l[3]):
		taxa[l[3]]['seqs'] += 1
		taxa[l[3]]['length'] += int(l[4])
	else: 
		taxa[l[3]] = {'seqs': 1, 'length': int(l[4])}
p.stdout.close()
print '  ' + str(len(recs)) + ' of all ' + str(allrecs) + ' sequences belonging to ' + str(len(taxa)) + ' taxa are retained.'


## Step 5: Choose a subset of representative taxa ##

if rank:

	# Subsample by group #
	
	print 'Subsampling...'
	for id in taxa:
		nowid = id
		while 1:
			if taxdumps[nowid]['rank'] == rank:
				taxa[id][rank] = nowid
				if groups.has_key(nowid): groups[nowid].append(id)
				else: groups[nowid] = [id]
				break
			else:
				if taxdumps[nowid]['parent'] == '1': break
				nowid = taxdumps[nowid]['parent']
				if not taxdumps.has_key(nowid): break
	
	group_nums = {}
	for id in groups:
		group_nums[id] = len(groups[id])
	
	rankppl = rank + 's'
	if rank[-7:] == 'species': rankppl = rank
	if rank[-5:] == 'genus': rankppl = rank[:-5] + 'genera'
	if rank[-6:] == 'family': rankppl = rank[:-6] + 'families'
	if rank[-5:] == 'class': rankppl = rank[:-5] + 'classes'
	if rank[-6:] == 'phylum': rankppl = rank[:-6] + 'phyla'
	
	print '  ' + str(len(group_nums)) + ' ' + rankppl + ' are contained in these taxa.'
	fout = open(rank + '_list.txt', 'w')
	for (id, num) in sorted(group_nums.items(), reverse=True, key=lambda x: x[1]):
		fout.write(id + '\t' + str(num) + '\t' + taxdumps[id]['name'] + '\n')
	fout.close()
	print '  A list of ' + rankppl + ' included in the BLAST database is saved as ' + rank + '_list.txt.'

	for gid in groups:
		nowsize = 0
		# Keep preferred taxons, if any.
		if lprefer:
			for id in groups[gid]:
				if lprefer.has_key(id):
					taxa[id]['in'] = 1
					nowsize += 1
					if nowsize >= size: break
		if nowsize >= size: continue
		# The 'representativity' is measured by the number of names a taxon has received from researchers.
		names = {}
		for id in groups[gid]:
			names[id] = taxdumps[id]['names']
		for (id, number) in sorted(names.items(), reverse=True, key=lambda x: x[1]):
			taxa[id]['in'] = 1
			nowsize += 1
			if nowsize >= size: break
	if size == 1: print '  Up to ' + str(size) + ' organism from each ' + rank + ' is sampled.'
	else: print '  Up to ' + str(size) + ' organisms from each ' + rank + ' are sampled.'

	# Keep certain taxonomic groups #

	if keep:
		keep_num = 0
		for id in taxa:
			nowid = id
			while 1:
				if nowid in keep:
					taxa[id]['in'] = 1
					keep_num += 1
					break
				else:
					if taxdumps[nowid]['parent'] == '1': break
					nowid = taxdumps[nowid]['parent']
					if not taxdumps.has_key(nowid): break
	
		keep_names = []
		for id in keep:
			if taxdumps.has_key(id):
				keep_names.append(taxdumps[id]['name'] + ' (' + id + ')')
		print '  ' + str(keep_num) + ' organisms from ' + ', '.join(keep_names) + ' are kept.'

else:

	# No subsampling #
	
	for id in taxa:
		taxa[id]['in'] = 1

# Print selected taxa and sequences #

(subtaxa, subrecs) = ({}, {})
for id in taxa:
	if taxa[id].has_key('in'):
		subtaxa[id] = taxdumps[id]['name']
for oid in recs:
	if subtaxa.has_key(recs[oid]['taxid']):
		subrecs[oid] = recs[oid]
if rank:
	print '  ' + str(len(subrecs)) + ' sequences from ' + str(len(subtaxa)) + ' organisms have been subsampled.'


fout = open('selected_organisms.txt', 'w')
for (id, name) in sorted(subtaxa.items(), key=lambda x: x[1]):
	srank = '0\tNA'
	if taxa[id].has_key(rank): srank = taxa[id][rank] + '\t' + taxdumps[taxa[id][rank]]['name']
	fout.write(id + '\t' + name + '\t' + srank + '\t' + str(taxa[id]['seqs']) + '\t' + str(taxa[id]['length']) + '\n')
fout.close()
print '  A list of subsampled organisms is saved as selected_organisms.txt.'

fout = open('selected_sequences.txt', 'w')
for oid in sorted(subrecs):
	fout.write(subrecs[oid]['gi'] + '\t' + subrecs[oid]['accn'] + '\t'  + subrecs[oid]['taxid'] + '\t' + subrecs[oid]['length'] + '\n')
fout.close()
print '  A list of subsampled sequences is saved as selected_sequences.txt.'


## Step 6: Create new database ##

print "SubDBer is ready to build a new database based on the subsampled sequences."
if not quiet: raw_input('Press any key to continue...')

fout = open('selected_gis.txt', 'w')
for oid in subrecs:
	fout.write(subrecs[oid]['gi'] + '\n')
fout.close()

if outfmt == 'gilist':
	print 'Creating GI list ' + outdb + '.gil...'
	p = subprocess.Popen('blastdb_aliastool -gi_file_in selected_gis.txt -gi_file_out ' + outdb + '.gil', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	if 'command not found' in out: sys.exit('Error: blastdb_aliastool not available. Check if you have ncbi-blast+ correctly installed.')
	print out.rstrip()
	os.remove('selected_gis.txt')	
	sys.exit('The task is completed.')

elif outfmt == 'alias':
	print 'Creating BLAST database alias...'
	p = subprocess.Popen('blastdb_aliastool -gilist selected_gis.txt -db ' + indb + ' -out ' + outdb + ' -dbtype ' + dbtype + ' -title ' + title, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	if 'command not found' in out: sys.exit('Error: blastdb_aliastool not available. Check if you have ncbi-blast+ correctly installed.')
	print out.rstrip()
	os.remove('selected_gis.txt')
	print 'Please keep ' + outdb + '.' + dbtype[0] + 'al and ' + outdb + '.' + dbtype[0] + '.gil at the same location or move to where the original BLAST database is located.'
	sys.exit('The task is completed.')

print 'Extracting sequences from the original database...'
p = subprocess.call('blastdbcmd -db ' + indb + ' -entry_batch selected_gis.txt > ' + outdb + '.fasta', shell=True)
os.remove('selected_gis.txt')

print 'Creating new ' + outfmt.upper() + ' database ' + outdb + '...'
if outfmt == 'blast':
	fout = open('selected_gi2taxids.txt', 'w')
	for oid in subrecs:
		fout.write(subrecs[oid]['gi'] + '\t' + subrecs[oid]['taxid'] + '\n')
	fout.close()
	p = subprocess.Popen('makeblastdb -in ' + outdb + '.fasta -parse_seqids -dbtype ' + dbtype + ' -out ' + outdb + ' -title ' + title + ' -taxid_map selected_gi2taxids.txt', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	for line in out.split('\n'):
		if not line: continue
		if line[0] == '\t': line = '  ' + line[1:]
		print '  ' + line
	p.stdout.close()
	os.remove('selected_gi2taxids.txt')
	os.remove(outdb + '.fasta')
	
print 'The task is completed.'

