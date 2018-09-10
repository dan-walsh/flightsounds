# flightsounds
## Overview
The **flight_event_extractor script**, which runs in the arcpy environment,  is designed to calculate the time and location of the closest approach of an aircraft to a sound station. It accomplishes this by linearly interpolating (in 3 dimensions) aircraft GPS points contained within input files provided by the user. The output is a csv file containing a list of events for a given station with associated location and velocity data, as well as flight metadata..

The csv file can then be used as an input to the **Flight_Event_to_SRCID script**, which builds a standard SRC_ID file using the output of the overflights script, which runs in the soundDB environment.

The resulting src_id file contains reference event numbers. When opened in SPLAT, each calculated aircraft approach for a given site will appear as a standardized rectangle centered on the event. The user must then manually splat each event , using the event id as a source number. The events can then be joined back to the overflights data using **SPLAT_append**, which joins the original csv flight event file to the results of splatting.

## Required Inputs and Outputs:

The input and output file names can be specified when running the scripts; names are suggested for consistancy.

### Input files:
    
#### allstations.shp
A file with zpoints describing the location of the sound stations (elevation in meters)
    	Must include following fields:
  	  	code - contains the four digit site code
   		 elevation - contains the elevation of the sound station in meters
    
#### flight_pts.shp
A file containing GPS data points describing the flight. Must include z-data in meters.
    	Must include the following fields:
   	 	altitude - contains the absolute MSL altitude of the plane at each point (int)
   	 ltime - contains the local time at which the point was recorded in the format 	"%Y/%m/%d %H:%M:%S"
   	 desc_ - contains a string with a unique flight identifier code
    
#### flight_lines.shp

A file containing lines generated by interpolating the points of each given flight. Must contain z data in meters. 
    Must include the following fields:
desc_ - the unique flight identifier, matched to the identifiers in the flight_pts file (string)


### Output Files
By default, output files are named SITE_overflights.csv, where SITE is the specific site code. They are saved directly in the selected working directory.

## Running the Scripts

### flight_event_extractor

*Note: This script **must run in a python 2.x environment on a computer running arcmap**. A yml file for building such an environment, named python3soundDBenv.yml, can be found in this repository. If the version of arcmap is different than Desktop 10.5, the following lines in the script need to be changed:*
```sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\arcpy")
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\bin")```

The ultimate goal of the first two scripts is to create two files: a srcID file containing predicted overflight events, and an overflights output file containing detailed information about the location and velocity of the plane at the time of closest approach.

To run the **flight_event_extractor** script, place the input shapefiles (allstations.shp, flight_pts.shp, and flight_lines.shp) in a single folder. The script will prompt the user to enter the path to this "working directory" when it is run. If the file names are changed the names entered in the script will need to be changed.

The script requires two additional inputs from the user: 

#### flight search radius 
The radius outward from the sound station within which flights will be assessed for point of closest approach (measured in meters). Values below 5000 are not reccomended as planes frequently travel these distances during the interval between points, meaning it is unlikely the script will have the two points it needs to calculate the closest approach.

#### station code 
The four digit code of the sound station you would like to calculate the closest approach to.

### Flight_Event_to_SRCID

*Note: This script **must run in a python 3.x environment**. A yml file for building such an environment, named python3soundDBenv.yml, can be found in this repository.*

This script takes the output of the flight_event_extractor, eliminates all closest approach events that do not have matching data in the sound database (accessed via soundDB), and saves the remaining events in SRCID format, allowing them to be manually quantified in splat. 
