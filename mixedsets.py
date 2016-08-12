# written by Sung-Hun Bae
import sys
import csv
from xml.dom.minidom import parse
from os import walk
from argparse import ArgumentParser

# merge chemical shifts less than a merge_threshold
def merge_peaks (sortedlist, merge_threshold) :
	i=0
	n=len(sortedlist)
	merged=[]
	neighbor=[ sortedlist[0] ]
	while i < n-1 :
   		if sortedlist[i+1] < sortedlist[i] + merge_threshold :
			neighbor.append ( sortedlist[i+1] )
		else:
			center= sum (neighbor) / len (neighbor)
			merged.append ( center )
			neighbor=[ sortedlist[i+1] ]
 		i += 1
	center= sum (neighbor) / len (neighbor)
	merged.append ( center )
	return merged

"""
	return min peak distance

	peaks: list of lists, [[v1,v2,..],[v1,v2,..],...]
	index: molecular index [v1,v2,..]
	overlap: peaks are considered overlapped if d is smaller than this
"""
def min_peak_distance (peaks,index,overlap) :
	# number of indices should be 2 or more
	if len(index) == 1 : 
  		return None
	# put a molecular index to each peak
	val=[]
	mol=[]
	for i in index :
		val.extend( peaks[i] )
		mol.extend( [i]*len( peaks[i]) )
	sortedIndex= sorted(range(len(val)), key=lambda x: val[x])
	val= sorted(val)
	mol= [ mol[x] for x in sortedIndex ]
	""" 
	remove overlapped peaks
	For a sorted peak list, val, measure distance, d,  
	between neighboring peaks. If d is larger than 
	overlap keep this distance otherwise replace 
	peaks' molecular index with -1
	"""
	prev= None
	for k, v in enumerate(val) :
		if prev is not None:
			d = v - prev
			if d < overlap : 
				mol[k] = -1
				mol[k-1] = -1
		prev = v

	# At least one peak must be present in the final set
	# so that it represents the molecule in the spectrum
	for i in index:
		if i not in mol:
			return None

	# return minimum peak distance
	mpd= []
	prev= None
	for k, v in enumerate(val) :
		if mol[k] != -1 :
			if prev is not None :
				d = v - prev
				mpd.append (d)
			prev = v
	return min(mpd)

def print_cluster ( number, l, cluster, d, label=None ) :
	if d :
	  	print "Set %2d (min peak distance: %6.3f ppm)" % (number,d), 
	else:
	  	print "Set %2d (min peak distance:        n/a)" % (number), 
	print len(cluster),"molecules"
	if not label :
		for i in sorted(cluster) :
			print "%4d" % i,
	else:
		for i in sorted(cluster) :
			print "%-20s" % label[i],
			for j, ppm in enumerate(sorted(l[i])) :
				print "%10.3f" % ppm,
				if j % 5 == 4 : print "\n%-20s" % ( label[i]),
			print
	print

# quality threshold clustering
def qtc (peaks,label,nmax,threshold,overlap):
	assigned= [ None ] * len(peaks) # reset cluster assignment
	cluster= [] # contains cluster member indices
	trial= [] # contains trial cluster memeber indices
	cluster_number= 1
	flag= True
	while assigned.count( None ) > 1 :
		cluster=[]
		for i, iidx in enumerate ( assigned ) :
			if iidx is None :
				trial= [ i ]
				flag = True
				while flag and assigned.count( None ) > len( trial ) :
					min_peak_d= None
					new_member= None
					for j, jidx in enumerate ( assigned ) :
						if ( jidx is None ) and ( j not in trial ) : 
							trial.append ( j )
							d = min_peak_distance ( peaks, trial, overlap )
							#print cluster,trial, d, min_peak_d, threshold, new_member
							if d and \
								(( min_peak_d is None ) or \
								( min_peak_d < d )):
								min_peak_d = d
								new_member = j
							trial.remove ( j )
					if new_member and min_peak_d > threshold and \
                       len(trial) < nmax :
						trial.append (new_member)
					else :
						flag = False

				# keep a set {A1,A2,..,A|G|} with maximum cardinality
				if len(cluster) < len(trial) :
					cluster = trial

		# assign a cluster number
		for i in cluster :
			assigned [i] = cluster_number;

		d = min_peak_distance (peaks, cluster, overlap)
		print_cluster (cluster_number, peaks, cluster, d, label )
		cluster_number += 1

	# while
  
	# give a set number for the remaining member
	for i, idx in enumerate( assigned ) :
		if idx is None :
			assigned [i] = cluster_number
			cluster= [i]
			print_cluster (cluster_number, peaks, cluster, None, label )
			cluster_number += 1




if __name__ == "__main__" :
	
	parser= ArgumentParser(description='NMR Chemical Shift Clustering')
	parser.add_argument('dirs', metavar='DIR', nargs='+',
 		help='Bruker NMR data directories')
	parser.add_argument('--merge',default=0.01,
		help='threshold to merge splitted peaks within a molecule (ppm)')
	parser.add_argument('--remove',default=0.02,
		help='threshold to remove overlapped peaks between molecules (ppm)')
	parser.add_argument('--cluster',default=0.3,
		help='threshold to cluster molecules (ppm)')
	parser.add_argument('--nmax',default=10,
		help='max number of molecules in a cluster')
	args= parser.parse_args()

	peaks= []
	label= []

	args_merge = float(args.merge)	
	args_remove= float(args.remove)
	args_cluster= float(args.cluster)
	args_nmax= int(args.nmax)

	print
	print "MixedSets by Sung-Hun Bae 3/18/2016"
	print
	print "Procedures and Parameters"
	print
	print "0. Given NMR data directories are walked through by mixedsets to search"
	print "   for a 'title' and a XML file for molecular ID and peak list,"
	print "   respectively. If the title file is empty, molecular ID uses"
	print "   the directory name instead." 
	print "1. A peak list defined in a XML file for each molecule"
	print "   is pre-processed in order to merge J splitted peaks."
	print "   [ threshold to merge   (--merge)  : %6.3f ppm]" % (args_merge)
	print "2. Peaks considered as overlapped between two molecules are"
	print "   removed from both peak lists during analysis of trial clusters."
	print "   At least one peak must remain for each molecule in a cluster."
   	print "   [ threshold to remove  (--remove) : %6.3f ppm]" % (args_remove)
	print "3. Molecules are clustered such that a distance between any two"
	print "   non-overlapped peaks be farther than a given threshold in ppm."
	print "   At least one unique and non-overlapped peak exists for each"
	print "   and every molecule in a cluster when mixed for screening."
	print "   [ threshold to cluster (--cluster): %6.3f ppm]" % (args_cluster)
	print "   [ max number of molecules in a cluster: %d]" % (args_nmax)
	print
	print "NMR Peak Data"

	for datadir in args.dirs :
		for (cur_dir, sub_dir, cur_files) in sorted(walk(datadir)) :
			title= None
			title_alt= None
			if "title" in cur_files :
				title_file= cur_dir+'/title'
				title_alt= '/'.join ( cur_dir.split('/')[-4:-2] )
				f= open(title_file,'r')
				title= f.readline().strip()
				f.close()

			if "peaklist.xml" in cur_files :
				peaklist_file= cur_dir+'/peaklist.xml'
				dom= parse (peaklist_file)
				pks= dom.getElementsByTagName('Peak1D')
				if not title:
					title= title_alt
				peaklist= []
				for p in pks :
					v= float( p.getAttribute('F1') )
					peaklist.append ( v )
				if peaklist :
					merged= merge_peaks ( sorted (peaklist), args_merge)
					print "%s" % ('/'.join(cur_dir.split('/')[-4:-2])),
					print "MolId", title,
					print " (%d peaks after merging %d)" % (len(merged), len(peaklist))
					for i, ppm in enumerate(sorted(merged)) :
						print "%10.3f" % ( ppm )
					print
					peaks.append ( merged )
					label.append ( title )
	print
	print "Quality Threshold Clustering"
	print
	qtc (peaks, label,args_nmax,args_cluster,args_remove)
