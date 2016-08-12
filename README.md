# mixedsets
Building a series of mixed sets of molecules based on their NMR chemical shifts.
Mixedsets use quality threshold clustering algorithm to cluster NMR chemical shifts dataset in order to resolve chemical shift overlap problem when you mix multiple compounds in a single NMR tube in a library screen.

# usage 
python mixedsets.py [-h] [--merge MERGE(default:0.01)] [--remove REMOVE(default:0.02)] [--cluster CLUSTER(default:0.3)] [--nmax NMAX(default:10)] DIR [DIR ...]

# procedures and parameters
1. Given NMR data directories are walked through by mixedsets to search
   for a 'title' and a XML file for molecular ID and peak list,
   respectively. If the title file is empty, molecular ID uses
   the directory name instead.
2. A peak list defined in a XML file for each molecule
   is pre-processed in order to merge J splitted peaks.
   [ threshold to merge   (--merge)  :  0.010 ppm]
3. Peaks considered as overlapped between two molecules are
   removed from both peak lists during analysis of trial clusters.
   At least one peak must remain for each molecule in a cluster.
   [ threshold to remove  (--remove) :  0.020 ppm]
4. Molecules are clustered such that a distance between any two
   non-overlapped peaks be farther than a given threshold in ppm.
   At least one unique and non-overlapped peak exists for each
   and every molecule in a cluster when mixed for screening.
   [ threshold to cluster (--cluster):  0.300 ppm]
   [ max number of molecules in a cluster: 10]

# example
```python mixedsets.py example/```

# example output
<pre>
NMR Peak Data
PL109A06/80 MolId PL109A06/80  (19 peaks after merging 21)
     0.703
     0.714
     0.725
     1.125
     1.136
     1.147
     1.158
     1.434
     1.445
     1.456
     2.500
     2.522
     2.533
     2.545
     7.652
     7.663
     7.678
     7.690
     8.235

PL109B06/181 MolId PL109B06/181  (5 peaks after merging 11)
     2.500
     6.529
     6.659
     7.114
     7.126

PL109C06/40 MolId PL109C06/40  (7 peaks after merging 9)
     2.500
     5.914
     6.807
     6.819
     7.027
     7.173
     7.187


Quality Threshold Clustering

Set  1 (min peak distance:  1.113 ppm) 2 molecules
PL109A06/80               0.703      0.714      0.725      1.125      1.136 
PL109A06/80               1.147      1.158      1.434      1.445      1.456 
PL109A06/80               2.500      2.522      2.533      2.545      7.652 
PL109A06/80               7.663      7.678      7.690      8.235
PL109C06/40               2.500      5.914      6.807      6.819      7.027 
PL109C06/40               7.173      7.187

Set  2 (min peak distance:        n/a) 1 molecules
PL109B06/181              2.500      6.529      6.659      7.114      7.126 
PL109B06/181
</pre>
