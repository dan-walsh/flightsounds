
# coding: utf-8

# In[1]:


import sys
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\arcpy")
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\bin")

import os
import datetime
import numpy as np
import pandas as pd
import csv
import io
import sys
import iyore
import soundDB


def makeSPLAT(filepath = r"V:\Noncanonical Data\2018 Concessionaire Overflight Analysis\2017 DENAUWBT Upper West Branch Toklat\DENAUWBT_overflights.csv",outputdir = r"C:\overflights\orig_nvspls"):

    # import the output of event extraction
    flight_events = pd.read_csv( filepath, sep=",")

        
   
    # make a list of unique values
    values = flight_events.StationID.unique()
    
    # creates a SRCID file for each unique StationID code in the inpuit document
    for value in values:
        # select all entries with a slingle code
        single_code = flight_events.loc[flight_events["StationID"] == value, "StationID": ]
        # create an output CSV file with just a subset of the entire dataset with matching codes
        single_code.to_csv(outputdir + value + ".csv", sep = ',')
        # create output filename
        out = outputdir + value + ".csv"
        # create output splat files using code_to_splat function below
        code_to_splat(out,value,outputdir)


def code_to_splat(file, code,outputdir):
    # import the output of event extraction
    flight_events = pd.read_csv( file,sep=",")


    # extract the times of closest approach for each event    
    flight_events["ClosestApch"] = pd.to_datetime(flight_events["ClosestApch"])

    print("Transcribing columns for " + code)

    # create an empty srcid dataframe
    srcid = pd.DataFrame()

    # column 1: nvsplDate
    flight_events["ClosestApch_date"] = flight_events["ClosestApch"].apply(lambda x: x.date())
    srcid['nvsplDate'] = flight_events["ClosestApch_date"]

    # column 2: hr
    flight_events["ClosestApch_time"] = flight_events["ClosestApch"].apply(lambda x: x.time())
    srcid['hr'] = flight_events["ClosestApch_time"].apply(lambda x: x.hour)

    # column 3: secs
    flight_events["ClosestApch_mins"]=flight_events["ClosestApch"].apply(lambda x: x.time())
    flight_events["ClosestApch_mins"]=flight_events["ClosestApch"].apply(lambda x: x.minute)

    flight_events["ClosestApch_secs"]=flight_events["ClosestApch"].apply(lambda x: x.time())
    flight_events["ClosestApch_secs"]=flight_events["ClosestApch"].apply(lambda x: x.second)

    # convert minutes to seconds
    flight_events["ClosestApch_mins"]=flight_events["ClosestApch_mins"].apply(lambda x: 60* x)

    # add seconds from above step to seconds
    srcid['secs']= flight_events["ClosestApch_mins"] +flight_events["ClosestApch_secs"]

    # column 4: len
    srcid['len'] = 60

    # column 5: srcID
    srcid['srcID'] = flight_events["EventID"].apply(lambda x: "1."+str(x).zfill(3))

    #column 6
    srcid['Hz_L']= 12.5

    #column 6
    srcid['Hz_U']= 2000

    #column 7
    srcid['MaxSPL'] = 0

    #column 8
    srcid['SEL'] = 0

    #column 9
    srcid['MaxSPLt'] = 0

    #column 10
    srcid['SELt'] = 0

    #column 11
    srcid['userName'] = "EXTRACTOR"

    #column 9
    srcid['tagDate'] = "1000-01-01 00:00:01"
    
    # column for sound data presence/absence -> will be populated after all processing
    srcid['data'] = 1
    
    # store initial frame length
    framelength = len(srcid)
    
    ###############
    ################
    # Check for audio data
    ###############
    print("Removing points without available audio for " + code)
    
    failure_count = 0
    
    ds = iyore.Dataset("Y:")
    
    # for every event in the srcid list
    for i in srcid.index:
        # extract time data
        ymd = srcid.iat[i,0]

        hour = srcid.iat[i,1]

        year = ymd.year
        month = ymd.month
        day = ymd.day
    
        #check availablity of a given event in the soundDB database
        available = [e for e in ds.nvspl(site=code, year=year, month=month, day=day, hour=hour)]

        # if nothing is available, set avaialable field to zero
        if(len(available)==0):
            srcid.ix[i,13]=0
            failure_count = failure_count+1
    
    
    

    # generate the file name
    filename = "SRCID_DENA" + str(flight_events.iloc[0,0])+ "eventExtractor.txt"

    print("Sorting and extracting data.")
    #sort the dataframe by date/time
    srcid = srcid.sort_values("nvsplDate")
    srcid = srcid.loc[srcid["data"] == 1]
    srcid = srcid.drop("data",axis=1)

    # export to csv with correct naming convention based on station code
    srcid.to_csv(outputdir + "SRCID_DENA" + code + ".txt",sep = "\t",index = False )

    # function for adding header lines to files
    def prepend(path, text):
        with open(path, 'r+') as f:
            body = f.read()
            f.seek(0)
            f.write(text + body);

    # add the splat format version header to the file
    prepend(outputdir + "SRCID_DENA" + code+".txt","%% SRCID file v20111005" + "\n")
    print("")
    print( str(framelength - failure_count) + " of "+ str(framelength) + " points had matching data for station id " + code +".")
        
filepath = input("Enter overflights calculator output file csv: ")


# if the user provided workspace is invalid, ask user again and again until valid path provided.
while not os.path.exists(filepath):
    print("")
    print(filepath)
    print("")
    print("Invalid path.")

    filepath = input("Enter valid path: ")

# print the new workspace value
print("")
print("File path changed to:")
print(filepath)
print("")

# check to see if it's a valid file by looking for EventID column. If it's not, ask for a new name.
flight_events = pd.read_csv( filepath, sep=",")

# check to see if file is valid. if not valid, ask for new file. check if new file exists. if it exists, check if valid. 
# if not exists, ask for new file. if it is valid, continue execution. 
while( 'EventID' not in flight_events.columns):
    print("Not a valid overflights output file.")
    filepath = input("Enter valid path: ")
    while not os.path.exists(filepath):
        print("")
        print(filepath)
        print("")
        print("Invalid path.")

        filepath = input("Enter valid path: ")
    flight_events = pd.read_csv( filepath, sep=",")

# set output directory by extracting directory from filepath using dirname
outputdir = os.path.dirname(filepath)+"/"
print("Output files will be saved to: " + outputdir)
print("")
input("Press enter to convert to SPLAT.")
print("")
# run scripts with validated input file and extracted output directory.
makeSPLAT(filepath , outputdir)
print("")
print("Splat conversion complete.")

