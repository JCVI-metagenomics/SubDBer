SubDBer
==========

Flexible taxonomy-based subsampling of sequence database for efficient and focused metagenomic analyses

Version: 0.1.1

Author: Qiyun Zhu (<qzhu@jcvi.org>)

Affiliation: J. Craig Venter Institute, La Jolla, CA, USA

Last updated: Oct. 4, 2015

License: [BSD 2-clause](http://opensource.org/licenses/BSD-2-Clause).

## Overview

Reference sequence databases are necessary for mapping, profiling, annotation and other metagenomic analyses. Larger databases provide higher accuracy at the cost of higher computational expense. We developed **SubDBer**, a simple and automated pipeline that allows researchers to build customized sequence databases that only include data of interest at reasonable comprehensiveness, based on resampling of large, standard BLAST databases. The program picks one or more representative taxa from each taxonomic group defined by the rank of users' choice. Users can further designate taxonomic groups to include, to exclude, and to escape subsampling (all taxa will be retained).

## Example

    python subDBer.py -in nt -out sub_nt -outfmt blast -within 2 -exclude 1117 -rank genus -size 1 -keep 816,838,1263

This command takes the NCBI nt database as input -> starts with all organisms from Bacteria (TaxID: 2) except for Cyanobacteria (1117) -> picks one representative organism per genus -> except for *Bacteroides* (816), *Prevotella* (838) and *Ruminococcus* (1263), in which all organisms are included -> creates a new BLAST database **sub_nt** that contains sequence data from the selected organisms.
