SpikeSort: semi-automated spike clustering software for matlab
Scott Grimes
Max Planck Institute - 2011

[Screenshot](SpikeSort/presentation/screenshot.png)

see presentation/presentation.pdf for more info

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
File List:

batchfilelist.txt
build_templates.m
Cluster.exe
cluster.m
cluster_consol.m
cluster_spc.m
combine_clusters.m
do_clustering.m
extract_data.m
find_temp.m
force_membership_wc.m
group_clusters.m
logfile.txt
logo.png
main_sorting.fig
main_sorting.m
mergefilelist.txt
nearest_neighbor.m
obtain_spikes.m
remove_cluster.m
run_cluster.m
test_ks.m
wave_features.m

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exporting from Spike2:

select save as: matlab format
select channels to export and timespans

the resulting .mat file should contain a structure for every channel

ex: testfile.mat
rat401_Ch1
rat401_Ch2
rat401_Ch3
rat401_Ch4

the following variables are required in order to run:
channel_structure.values
channel_structure.length
channel_structure.start
channel_structure.interval

you can cluster any number of channels together in the same file

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Single File:
1) select the .mat file you wish to cluster 

2) the box underneath the filename is the standard 
   deviation amplitude threshold for spike detection
   (default value is 5mu)

3) choose the length of time in ms you wish to obtain
   from the moment the amplitude threshold is crossed.
   left is .3ms and right is .7ms by default

4) select the meathod of clustering (SPC or PCA)
   if you know the number of clusters ahead of time you can
   input the value to the right of the selection box,
   otherwise leave the value at 0 for auto-selection

5) hit "go" undereath the filename

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Batch Files:
1) open "batchfilelist.txt" and write each file name including the
   ".mat" extention on a seperate line.

2) the box underneath the filename is the standard 
   deviation amplitude threshold for spike detection
   (default value is 5mu)

3) choose the length of time in ms you wish to obtain
   from the moment the amplitude threshold is crossed.
   left is .3ms and right is .7ms by default

4) select the meathod of clustering (SPC or PCA)
   if you know the number of clusters ahead of time you can
   input the value to the right of the selection box,
   otherwise leave the value at 0 for auto-selection

5) hit "batch clustering"

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Plotting:
Type the numbers of the clusters you wish to plot in the box seperated
by a space. 0 is default, and plots all clusters. the checkbox is checked
by default. unchecking it will plot every spike individually, leaving it checked
will only plot the cluster averages and standard deviation.

Raw data can be plotted in the bottom left corner.

The spike timelines can be plotted in the bottom left corner. 0 is default, and plots all clusters
on one graph. You can select individual clusters or multiple clusters to be plotted by typing the
cluster number into the box on the right, with each number seperated by a space.

the inter-spike interval histograms can be plotted in the bottom left corner.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Saving/Loading:
Type in the name you wish to save/load the data as, along with a ".mat" extention. 

saved data contains 4 variables:
ts time vector from start to end, spaced with interval
spikes a NxP vector of N spikes of P length

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Merging Files:
1) open mergefilelist.txt and type each filename to be merged on a seperate line (with ".mat" extention)

2) choose the final filename in the gui and click merge files

3) the files from previous runs will be combined into one, and can be clustered together

*note- only files which have already run through the spike detection portion can be merged. (batch files also accepted)








~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exporting to Spike2:
1) after loading a clustered file, type the desired name into the save_file box
2) click export to spike2
3) the .txt file is saved, and is loaded using the "importation.s2s" script in spike2
4) verify that the variable x in "num[x]" in the script is the same length as your spikes
5) run the script and select your exported .txt file









