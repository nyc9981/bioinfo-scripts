#!/usr/bin/env python
# modified from /home/pskene/bin/py_peak_calling.py


import click

@click.command()
@click.argument('bedgraph', type=click.File('r'))
@click.option('--threshold', default=1, help='minimum coverage to be included in peaks')
@click.option('--min_length', default=20, help='minimum length of peaks')
@click.option('--inter_peak_distance', default=1, help='Inter-peak distance cutoff')
@click.option('--merge_close_peaks', default=True, help='Merge close peaks')
@click.option('--keep_highest_close_peak', default=False, help='Keep the highest close peak')
@click.option('--max_length', default=10000, help='maximum length of peaks')
@click.option('--generate_ID', default=True, help='Automatic generation of an integer list of peak IDs')
@click.option('--output_name', default=None, help='User defined peak file name')
@click.option('--delete_overlap_bed', default=None, help='bedfile whereby any peaks that overlaps with this bedfile will be discarded')
def py_peak_calling(bedgraph, threshold, min_length, inter_peak_distance, merge_close_peaks, keep_highest_close_peak, max_length, generate_ID, output_name, delete_overlap_bed):
    """
    - need to install a more up-to-date varsion of bedtools before invoking Jupyter
      type: module load bedtools/2.21.0
	(1) filters bedgraph based on threshold;
	
	(2) merges adjacent basepairs that are over threshold;
	
    (3) retains peaks that satisfy min/max length criteria; 
	
	(4) merges any peaks that are closer than the inter-peak distance cutoff -or- alternatively keeps just the highest peak (this is beta functionality)
	
    - max length is typically defaulted to be very large
    - outputs a bed file (default col4 is the sum of the bedgraph scores; sorted by chrom;start;stop)
    - generate ID: will auto generate a integer list as a ID number (1... number of peaks). This will be reported as column 4 and the bedgraph scores will be shifted to column 5 as per standard bed format
    - note the peak score for merged peak is the *just* the sum of the two individual peaks not the 
    total score in the merged region (i.e. there could be some sub-threshold scores in the intervening space that won't be included)
    -assumes bedgraph in standard format <chr> <start> <stop> <score>
    -output_name = option for user defined name (type with '...'), otherwise will generate name bedgraph_peaks.bed
    -delete_overlap_bed = option to add path to bedfile (as string), whereby any peaks that overlap this bed file will    be discarded
    """
    
    import pybedtools
    import glob
    from pybedtools import BedTool
    import pandas as pd
    import csv
    
    if merge_close_peaks==keep_highest_close_peak:
        return 'Exiting... merge_close_peaks and keep_highest_close_peak set the same'
    
    #generate name for output
    bedgraph_name = glob.glob(bedgraph)
    
    if output_name:
        filename = output_name
        
    else:
        filename = bedgraph_name[0].replace('.bg', '_peaks.bed')
        
    print 'input bedgraph file: ' + bedgraph_name[0]
    print 'output filename: ' + filename
    
    #import data as BedTool
    data = BedTool(bedgraph) 
    
    #retains intervals above threshold
    above_thresh = data.filter(lambda b: float(b.name) >= threshold) 
    
    #merge adjacent above threshold regions and sum bedgraph scores (assumes bedgraph score in col 4)
    #by increasing d value can allow for 
    merge_regions= above_thresh.merge(d=0, c=4, o='sum' )
    
    #filter based on length criteria
    peaks = BedTool(merge_regions.filter(lambda x: len(x) >= min_length and len(x) <= max_length))
    
    #print 'number of regions identified before merging or filtering: ' + str(peaks.count())
    
    if merge_close_peaks==True:
        #merge the bonafide peaks if they they are shorter than the inter peak distance and sum scores and sort
        print 'merging peaks that are closer than: ' + str(inter_peak_distance)
        merge_peaks = peaks.merge(d=inter_peak_distance, c= 4, o='sum').sort()
        
    if keep_highest_close_peak==True:
        #need to read each line to find close peaks and throw away the one with the lowest score out of the two
        print 'entering loop'
        
        peaks.saveas('temp_input.bed')
        
        print 'before keeping highest, number of regions identified: ' + str(BedTool('temp_input.bed').count())
        
        last_line = [str(item) for item in (BedTool('temp_input.bed').to_dataframe().tail(n=1).iloc[0,:].tolist())]

        with open('temp_input.bed') as myfile:    
            with open('test_output.bed', 'w') as output:
                file_output = csv.writer(output, delimiter='\t')

                prev_line = None

                for line in csv.reader(myfile, delimiter='\t'):

                    if prev_line is None:
                        prev_line = line
                    elif float(prev_line[2])+float(inter_peak_distance) <= float(line[1]):
                        file_output.writerow(prev_line)
                        prev_line = line
                    else:  
                        if float(prev_line[3]) < float(line[3]):
                            prev_line = line
                        # modified by Xin
                        else:
                            file_output.writerow(prev_line)
                            prev_line = line

                if line==last_line:
                    file_output.writerow(prev_line)
                    
            merge_peaks = BedTool('test_output.bed')
   
    print 'number of peaks found: ' + str(merge_peaks.count())
    
    if delete_overlap_bed!=None:
        print 'delete_overlap_bed provided: ' + delete_overlap_bed
        merge_peaks = merge_peaks.intersect(b=delete_overlap_bed, v=True)
        print 'number of peaks retained: ' + str(merge_peaks.count())
    
    if not generate_ID:
        print 'saving sorted peak bed file with no ID'      
        merge_peaks.saveas(filename)        
    else:
        print 'saving sorted peak bed file with ID names'
        
        #change to pandas dataframe
        DF_peaks = merge_peaks.to_dataframe()
        
        #insert new column with id: 1.... # of peaks
        DF_peaks.insert(3, 'id', ['id' + str(item) for item in range(1, (len(DF_peaks)+1))])

        #save output
        DF_peaks.to_csv(filename, sep = '\t', header = False, index = False)
        
    return 'Finished'
    
if __name__ == '__main__':
    py_peak_calling()