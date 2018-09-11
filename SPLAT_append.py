
# coding: utf-8

# In[6]:


import sys
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\arcpy")
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\bin")

import os
import datetime
import numpy as np
import pandas as pd
import soundDB
import iyore
ds = iyore.Dataset("K:")
import matplotlib.pyplot as plt

# SPLAT_append combines a srcID file output by Flight_Event_to_SRCID  with the output of the flight_event_extractor
# script. The srcID file should contain both the output of the Flight_Event_to_SRCID script and the events manually
# annotated in SPLAT. 
#
#
# The resulting csv file contains a list of all overflight events with associated flight data, such as flight number, elevation,
# climb/descent rate, ground speed, etc. 
#
# example function call: SPLAT_append("C:/overflights","UWBT_overflights.csv","SRCID_DENAUWBT.txt")
#
# 




# change display options to make all rows and columns visible for debugging
pd.options.display.max_columns = 999
pd.options.display.max_rows = 999


def SPLAT_append(workingdirectory, flighteventfile, splatfile):
    # import the output of event extraction (flight_event_extractor script)
    flight_events_calculated = pd.read_csv( workingdirectory + "/" + flighteventfile ,sep=",")
    
    # extract the station code for use in output file name
    code = flight_events_calculated.StationID[0]
    
    # import srcid file (see description above)
    src = pd.read_csv( workingdirectory + "/" + splatfile,sep="\t", header = 1 )

    # separate extractor (calculated using script) vs  events manually splatted using AMT into two dataframes
    rows = src['userName'] == 'EXTRACTOR'
    # events calculated by script
    src_c = src[rows]

    # events splatted by user
    rows = src['userName'] != 'EXTRACTOR'
    src_s = src[rows]

    # append audibility column to user splatted events. When splatting, the pre-period digit 1 or 0 is used to indicate 
    # whether an event is visible in the spectrogram or not visible, respectively. 
    
    # make copy to avoid error message "copyfromslice"
    src_s = src_s.copy()
    src_s["audible"] = src_s["srcID"] > 1


    # a function to standardize the id numbers by removing the audibility prefix.
    def simpleeventnumber(item):
        item = round(item,3)
        if( item < 1):
            item = item + 1

        item = str(item)    
        item = item.split('.', 1)[1]

        return item
    # apply the function to the splatted and calculated events (since all calculated events will be 1.xxx but some inaudible splatted
    # will be 0.xxx)
    src_s = src_s.copy()    
    src_s['trueID'] = src_s.srcID.apply(simpleeventnumber)
    src_c = src_c.copy()    
    src_c['trueID'] = src_c.srcID.apply(simpleeventnumber)

    # set both indexes to the same column

    src_s.set_index('trueID')

    src_c.set_index('trueID')

    # and join the two tables on the trueID column, producing the final output csv.
    splat_correct = src_c.merge(src_s,how='inner',on = 'trueID')

    # number of non-matching events is the difference in length between splatted data and the calculated data
    print("Tables joined. Number of non-matching event IDs: " + str(abs(splat_correct.shape[0]-src_s.shape[0])) + "\n")

    # output only the matching events, removing events that were splatted with an incorrect number.
    output = src_s[src_s['trueID'].isin(splat_correct.trueID)]
    # chang true id to integer rather than a string
    output['trueID'] = src_s.trueID.apply(int)
    
    # count how many events there are and print the result for user output
    count = (str(len(output)))
    print(count + " matching events exported.")
    print("")

    # join "output" (the error-checked srcid file) to the output of the flight event extractor
    joined = output.merge(flight_events_calculated,how='inner',left_on = 'trueID',right_on = 'EventID')

    # and export the result as a csv
    joined.to_csv(workingdirectory + "/DENA" + code + ".txt")
    
    print("File exported to " + workingdirectory + "/DENA" + code + "_joined.txt")
    #######################################################
    #######################################################
    

# these three blocks allow user to change names of input files and working directory
workingdir = "c:/overflights"
print("Current working directory: " + workingdir)
change = input("Change working directory?")
if(change == "y"):
    workingdir = input("Enter new workspace directory: ")

    # if the user provided workspace is invalid, ask user again and again until valid path provided.
    while not os.path.exists(workingdir):
        print("Invalid path.")
        print("")
        workingdir = input("Enter new workspace directory: ")
    # print the new workspace value
    print("Workspace changed to:")
    print(workingdir)
    print("")

    
splatfile = input("Enter name of srcid file: ")
# if the user provided workspace is invalid, ask user again and again until valid path provided.
while not os.path.exists(workingdir +"/"+ splatfile):
    print("Invalid path.")
    print("")
    splatfile = input("Enter new srcid file name: ")
    print("")


flight_event_file = input("Enter name of flight event csv file: ")
# if the user provided workspace is invalid, ask user again and again until valid path provided.
while not os.path.exists(workingdir +"/" +flight_event_file):
    print("Invalid path.")
    print("")
    flight_event_file = input("Enter name of flight event csv file: ")
    print("")

# run the script with the provided inputs
SPLAT_append(workingdir,flight_event_file,splatfile)

