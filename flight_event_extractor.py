
# coding: utf-8

# In[6]:




import sys
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\arcpy")
sys.path.append(r"C:\Program Files (x86)\ArcGIS\Desktop10.5\bin")



import arcpy
import os
import datetime
import numpy as np



# Find point of closest approach.
# 
# Input: flight track shapefile (a line), sound station single point location with z data
# 
# output: a point containing the closest point of approach of the line to the sound station (no z data)

def closest_approach(station_loc, flight_ln, flight_path):

    
    # rename variables for local use
    sound_station_loc = station_loc #
    flight_line = flight_ln #
    flight_pts = flight_path #
    

    ##################################################
    # near 3d used to append distance to inputted sound station to data table for each flight track point
    #
    # copy data table to prevent changes to original file
    arcpy.CopyFeatures_management(flight_pts, "in_memory//flight_pts")
    
    # check out the exenstion
    arcpy.CheckOutExtension("3D")
     # apply the extension
    arcpy.Near3D_3d("in_memory//flight_pts",sound_station_loc,location="LOCATION")
    # check in the extension
    arcpy.CheckInExtension("3D")
    ###################################################


    #########################################################
    # append nearest approach of line to sound station to the line shapefile
    #
    # first, copy features to prevent overwriting
    arcpy.CopyFeatures_management(flight_line, "in_memory//flight_line")
   

    
    
    # then use near3d

    # check out the exenstion
    arcpy.CheckOutExtension("3D")
     # apply the extension
    arcpy.Near3D_3d("in_memory//flight_line",sound_station_loc,location="LOCATION")

    ########################################################


    ########################################################
    # create point of closest approach
    ###########################################################
    # extract the x and y coordinates from "in_memory//flight_line"
    # and use them to create a new point which will be passed to the extract_by_proximity function
    xcoordinate = []
    with arcpy.da.SearchCursor("in_memory//flight_line", ['NEAR_FROMX']) as cursor:
        for row in cursor:
            xcoordinate.append(row[0])
        
    ycoordinate = []
    with arcpy.da.SearchCursor("in_memory//flight_line", ['NEAR_FROMY']) as cursor:
        for row in cursor:
            ycoordinate.append(row[0])

    zcoordinate = [0]

    # create a point object, make it into a geometry, then add it to the shapefile
    ptList =[[xcoordinate[0],ycoordinate[0],zcoordinate[0]]]
    pt = arcpy.Point()
    ptGeoms = []
    for p in ptList:
        pt.X = p[0]
        pt.Y = p[1]
        pt.Z = 0
        ptGeoms.append(arcpy.PointGeometry(pt))

    # export the new point to shp, define the coordinate system to the same as flight line
    arcpy.CopyFeatures_management(ptGeoms, "in_memory//closest_pt")
    arcpy.DefineProjection_management ("in_memory//closest_pt", "in_memory//flight_pts")
    
    # add an altitude field to closest_pt
    arcpy.AddField_management("in_memory//closest_pt","alt", "SHORT")
    

    # convert into a 3d point feature for use by near3d
    arcpy.FeatureTo3DByAttribute_3d("in_memory//closest_pt", "in_memory//closest_pt_3d",'alt')
    closest_point = "in_memory//closest_pt_3d"
    


    # check in the extension
    arcpy.CheckInExtension("3D")

    # assign and return output of function
    output_data = ["in_memory//flight_pts",closest_point]
    return output_data;



# extract_by_proximity_pts - this function selects the two closest GPS points to the point of closest approach of a flight
# 
#    Input files:    .shp with z data point cloud 
#                    .shp containing a single point
#                    
#    Output files:   .shp containing n points from the point cloud
#    
#    Parameters:     n    number of points to consider
#                     
# Extracts the n closest points from the point cloud to the single point and saves them as a new shapefile, retaining all
# original fields.
# 



def extract_by_proximity(multipt, targetpt, npts, output_file):
    #set local vars
    n=npts

    # copy in_memory files to local variables for the function
    arcpy.CopyFeatures_management(multipt, "in_memory//test_cloud_copy")
    arcpy.CopyFeatures_management(targetpt, "in_memory//singlepointcopy")
    single_pt_copy = "in_memory//singlepointcopy"
    test_cloud_copy = "in_memory//test_cloud_copy"
    


    ##########################################################
    # use near 3d to append distance data to the pt_cloud points

    # check out the exenstion
    arcpy.CheckOutExtension("3D")
    # apply the extension
    arcpy.Near3D_3d(test_cloud_copy,single_pt_copy, angle="ANGLE")
    # check in the extension
    arcpy.CheckInExtension("3D")
    

                                                         
    ##########################################################
    # sort the resulting output shp table by distance
    arcpy.Sort_management(test_cloud_copy,"in_memory//test_cloud_copy_sorted",sort_field="NEAR_DIST3 ASCENDING", spatial_sort_method="UR")
    test_cloud_copy_sorted = "in_memory//test_cloud_copy_sorted"
    
                                  
      
    
    #####################################################################
    # extract the lowest n table entries and save as the final output file
    arcpy.Select_analysis(test_cloud_copy_sorted,output_file, '"FID"<' + str(n+1))
    #print("Proximity analysis: Complete. Closest " + str(n) + " points extracted.")
    return output_file;





# 3pt_field_interpolation - interpolates GPS data from two closest points to point of closest approach to estimate speed, heading.
# and other metrics
# 
#    Input files:    .shp with z data containing two points and the distance of each point from the middle point
#                    .shp containing a single point on a line between the two existing points
#                    
#    Output files:   .shp containing the single point and several new attribute table fields (always elevation, others optional)
#    
#    Parameters:     boolean "calc_speed"
#    
#   
#                     
# Interpolates altitude of intermediate point using a weighted average (linear) interpolation.

# In[4]:


def flight_interpolation(twopts, singlepoint, sound_station_loc, counter):

    # initialize note variable
    note = ""
    # name the needed shapefiles
    two_pts = twopts
    single_point = singlepoint
    
    # uncomment to outputthe  two selected points for troubleshooting
    #arcpy.Select_analysis(two_pts,"twopoints.shp")

    # read the distances from the two points and save to list 'point_distances'
    point_distances = []
    with arcpy.da.SearchCursor(two_pts, ['NEAR_DIST']) as cursor:
        for row in cursor:
            point_distances.append(row[0])



    # read the altitudes from the two points and save to list altitudes  
    altitudes = []
    with arcpy.da.SearchCursor(two_pts, ['altitude']) as cursor:
        for row in cursor:
            altitudes.append(row[0])        

      # read the times from the two points and save to list 'altitudes'  
    times = []
    with arcpy.da.SearchCursor(two_pts, ['ltime']) as cursor:
        for row in cursor:
            times.append(datetime.datetime.strptime(row[0], "%Y/%m/%d %H:%M:%S" ))            
        
    # calculate the total distance between the points
    total_dist = round(point_distances[0] + point_distances[1],0)

    ######################################
    # ALTITUDE CALCULATION - calculates the altitude at the closest approach point
    ######################################

    # calculate the total elevation change
    alt_change = altitudes[1]-altitudes[0]
    
    # eliminate points with total_dist of zero to prevent division by zero; add note if changed
    if(total_dist == 0):
        total_dist = 100000
        note = note + "Zero distance - duplicate points; "        
    
    # calculate the altitude change between point 0 and the point of interest
    alt_change_to_pt = alt_change*(point_distances[0]/float(total_dist))
    # calculate the final altitude of the point of interest by adding to the base altitude    
    point_alt = alt_change_to_pt + altitudes[0]
    
    # save altitude of closest point to pass to results
    altatclosestapchm = point_alt
    



    #######################################
    # CLOSEST APPROACH TIME - determine the time of the closest approach to the sound station
    #######################################

    #total time change
    time_change = (times[1]-times[0]).total_seconds()
    
    # add error message for zero time change, change to high number to ensure removeal by "unrealistic speed" test below
    if(time_change == 0):
        time_change = 10000000
        note = note + "Zero time change; "     

    #calculate the time change between point 0 and the point of interest
    time_change_to_pt = datetime.timedelta(seconds = time_change*((point_distances[0]/float(total_dist))))

    # calculate the actual closest approach time
    closest_apch_time = str(time_change_to_pt + times[0])
    
    # remove microseconds from approach time
    closest_apch_time = closest_apch_time.split('.', 1)[0]


    #########################################
    # average speed between points - average speed of the plane between the two closest point
    ##################################
    avg_speed = abs(round(total_dist / time_change,0))

  
    
    # measure the closest approach distance
    arcpy.PointDistance_analysis(single_point, sound_station_loc, "in_memory//distance")
    
    # read the 2d distance value from the table AND SAVE TO VARIABLE "DISTANCE"
    with arcpy.da.SearchCursor("in_memory//distance", ['distance']) as cursor:
        for row in cursor:
            distance = row[0]    
 
        
    # read the altitude of the station from the input file
    with arcpy.da.SearchCursor(sound_station_loc, ['elevation']) as cursor:
        for row in cursor:
            elevation = row[0]    
     
    
    # subtract the two elevations and use trig to calculate actual distance
    ele_diff = point_alt - elevation
    closest_dist3d = round(np.sqrt((distance*distance)+(ele_diff*ele_diff)),0)
        
        
    # extract the tail number and flight number    
    flight_id = []
    with arcpy.da.SearchCursor(two_pts, ['desc_']) as cursor:
        for row in cursor:
            flight_id.append(row[0])    
    # extract the station ID
    station_id = []
    with arcpy.da.SearchCursor(sound_station_loc, ['code']) as cursor:
        for row in cursor:
            station_id.append(row[0])    

    # determine climb or descent
    climbrate = round(-1*(alt_change_to_pt / int(time_change)),2)

    ##################################################
    # measure angle of flight, convert to compass bearing
    ##################################################
    # first, read the angle value from the twopts file:
    with arcpy.da.SearchCursor(twopts, ['NEAR_ANG_H']) as cursor:
        for row in cursor:
            angle = row[0] 
    

    # ensure that angle is on positive side
    if(angle <0):
        angle = angle + 180
    # convert geometric angle to north angle    
    if(angle <= 90):
        angle = 90-angle
        
    if(angle > 90):
        angle = 180 - angle


    # extract x coordinates of the two points
    xcoord = []
    with arcpy.da.SearchCursor(two_pts, ['SHAPE@X']) as cursor:
        for row in cursor:
            xcoord.append(row[0])    

    # orient the vector in the direction of motion based on the time data
    if(time_change > 0):
        if(xcoord[0]>xcoord[1]):
            angle = angle + 180
    if(time_change <= 0):
        if(xcoord[0]<xcoord[1]):
            angle = angle + 180   

    # round the final angle value
    angle = round(angle,1)
    
    #########################################################
    # extract x and y coordinates of closest approach point
    ########################################################

    # set the output CS (coordinate system) to ensure data will plot in arcmap (alaska albers 2011)
    outCS = arcpy.SpatialReference(102117)    
    arcpy.Project_management(single_point, "out.shp",outCS)

    # extract all x coordinates
    xcoordca = []
    with arcpy.da.SearchCursor("out.shp", ['SHAPE@X']) as cursor:
        for row in cursor:
            xcoordca.append(row[0]) 
    # extract all y coordinates
    ycoordca = []
    with arcpy.da.SearchCursor("out.shp", ['SHAPE@Y']) as cursor:
        for row in cursor:
            ycoordca.append(row[0]) 

    # print the spatial referenc of the output for troubleshooting
    desc = arcpy.Describe("out.shp")
    spatialRef = desc.spatialReference
 
    # Print the spatial reference name
    print spatialRef.Name

    ################################################
    # Function output:
    ################################################
    
    # create a list with all output variables
    output_list = [str(station_id[0]),str(flight_id[0]),str(counter).zfill(3),str(closest_dist3d),str(closest_apch_time),str(avg_speed*3.6),str(np.abs(time_change)),str(total_dist),str(climbrate),str(angle),note,str(xcoordca[0]),str(ycoordca[0]),altatclosestapchm]

    # print a summary of the event
    print("")
    print("Event Summary:")
    print("")
    print("Flight ID                       :  " + str(flight_id[0]))
    print("Station ID                      :  " + str(station_id[0] ))
    print("Distance at Closest Approach (m):  " + str(closest_dist3d))
    print("Time of Closest Approach        :  " + str(closest_apch_time))
    print("Average interval speed (km/hr)  :  " + str(avg_speed*3.6))
    print("Time Between Pts (s)            :  " + str(np.abs(time_change)))
    print("Distance Between Pts (m)        :  " + str(total_dist))  
    print("Climb Rate (m/s)                :  " + str(climbrate))  
    print("Vector (deg)                    :  " + str(angle))  
  
    # warn if speed is less than 50 km/hr
    print("")
    if np.abs(avg_speed*3.6) < 50:
        print("Unrealistic speed, please review input tracks.")
    # add any error codes accomulated while running the function
    if( note != ""):
        print("Error codes: " + note)
    # or not
    if( note == ""): 
        print("No data validation errors.")
    return output_list;


# single extraction uses the location of a specific station, GPS point cloud, and a single flight line to perform the interpolation for 
# this function combines the upper three functions into one function, which is called by the extract_events function
def single_extraction(station_loc, flight_ln, flight_pts,counter):
    # first, test to make sure there are two points in flight_pts
    # if not, return a null output and print the error message.

    global ouput
    x = int(arcpy.GetCount_management(flight_pts).getOutput(0))
    if( x < 2):
        print("Extraction Failed: < 2 points within search radius.")

        output=[0]
     ##############################################################
    # applies the three main functions to calculate the flight attributes,
    # sends the final results to the extract_events function to be written to output file
    else:
        cpout = closest_approach(station_loc, flight_ln, flight_pts)
        flight_pts = cpout[0]
        closest_pt = cpout[1]
        out = extract_by_proximity(flight_pts,closest_pt,2,"in_memory//outfile")
        
        # run flight interpolation to produce the final function output passed to extractevents function
        output = flight_interpolation(out, closest_pt,station_loc, counter)
        
    return output


# Input file requirements:
#     
# allstations
# 
# A file with zpoints describing the location of the sound stations (elevation in meters)
#     Must include following fields:
#     
#     code - contains the four digit site code
#     elevation - contains the elevation of the sound station in meters
#           all 3s station files used should contain z data (z-point) describing the elevation of each station
#            allstations = "3d_stations.shp"
#     
# Flight_pts
#     
# A file containing GPS data points describing the flight. Must include z-data in meters.
#     Must include the following fields:
#        
#     altitude - contains the absolute MSL altitude of the plane at each point (int)
#     ltime - contains the local time at which the point was recorded in the format "%Y/%m/%d %H:%M:%S"
#     desc_ - contains a string with a unique flight identifier code
#     
# Flight_lines
# 
# A file containing lines generated by interpolating the points of each given flight. Must contain z data in meters. 
#     Must include the following fields:
#     
#     desc_ - the unique flight identifier, matched to the identifiers in the flight_pts file (string)


    

# extract_events: Where it all comes together.
#
# extract_events will take as inputs flight points and flight tracks. It will iterate through all the tracks
# (as defined by the desc_ field) and run single_extraction on all of them. Single extraction will return a 
# single line of a csv file, which will be appended by extract_events onto the output csv file.
#

def extract_events(code,bufferdist,flight_pts = "2017trackpts.shp",flight_lines = "2017tracklines.shp", outputdir = "C:\\overflights\\", allstations = "3d_stations.shp"):

    # import needed libraries
    import csv 
    import arcpy
    
    # hardcoded references below replaced by parameterization of extract_events
    #outputdir = "C:\\overflights\\"
    #flight_pts = "2017trackpts.shp"
    #flight_lines = "2017tracklines.shp"
    #allstations = "station_loc3d.shp"
    

    # test to make sure station location data contains z points. If not, provide warning message.
    desc = arcpy.Describe(allstations)
    if str(desc.hasZ) == "False":
        print("Sound station location file does not contain z data. Please use \"Feature to 3d by attribute \" tool to add z data in meters.")
        return
    arcpy.Select_analysis(allstations,"in_memory//station",'"code"= \'' + str(code) + '\'')
    station = "in_memory//station"
    
    # check to make sure station exists
    if(int(arcpy.GetCount_management(station).getOutput(0)) != 1):
        print("No matching stations. Please check input code. \n")
        return

    # takes the user bufferdist input and appends meters to it, required for input to buffer analysis arcpy function
    # bufferdist should be set to somewhat larger than the desired search distance, especially for small values; this
    # is because a track may come much closer than the closest two points
    bufferdist = str(bufferdist) + " meters"

    # perform the buffer and use results to clip flight points and flight lines
    arcpy.Buffer_analysis(station, "in_memory//stationbuffer",bufferdist)
    arcpy.Clip_analysis(flight_pts, "in_memory//stationbuffer", "in_memory//clippedtrackpts")
    arcpy.Clip_analysis(flight_lines, "in_memory//stationbuffer", "in_memory//clippedtracklines")



    # function to extract the unique flight id values from the table (field: desc_) and save to variable flight_ids
    # this list will be used to iterate through all flights in the input table
    def unique_values(table, field):
        with arcpy.da.SearchCursor(table, [field]) as cursor:
            return sorted({row[0] for row in cursor})
    # run the function above to get the list
    flight_ids = unique_values("in_memory//clippedtracklines", "desc_")
    
    # output the num of pts that will be processed
    print(str(len(flight_ids))+ " points will be processed.")
    print("")

    # initialize the output csv file
    with open( outputdir + "/" + code + '_overflights.csv','wb') as textfile:
        writer = csv.writer(textfile, delimiter=',')
        writer.writerow(["StationID","FlightID","EventID","ClosestDist3d","ClosestApch","AveSpeedKPH","TimeChange","TotalDist","ClimbRate","Vector","Note","xcoord","ycoord","altmclosedist"])
    
    # following process carried out for every track. x used to keep track of how many have been processed
        global x
        x=1
        # for each unique flight value
        for value in flight_ids:

            # define the sql expression used to select the desired fields
            expression = '"desc_"= \'' + str(value) + '\''
            
            # print the name of the file currently being processed
            print("--------------------------------------------------------------- \n")
            print( str(x) + "/"+ str(len(flight_ids)) + "                            Processing " + str(expression))
            
            
            # select the an individual track and its lines by id, 
            arcpy.Select_analysis("in_memory//clippedtrackpts","in_memory//singletrackpt",expression)
            arcpy.Select_analysis("in_memory//clippedtracklines","in_memory//singletracklines",expression)
            # peform the extraction for the selected event
            outputdata = single_extraction(station, "in_memory//singletracklines", "in_memory//singletrackpt",x)

            # increment the counter
            x = x+1


            # only write the output to the csv is it isn't equal to zero (which only happens if the interpolation script detects unreasonable outputs)
            if(outputdata[0] != 0):
                writer.writerow(outputdata)
    # close text file and clean up memory
    textfile.close()
    arcpy.Delete_management("in_memory")
    print("In_memory is cleared")

    return

##############################################################
##############################################################
##############################################################
### USER INTERFACE AND SCRIPT EXECUTION BELOW
##############################################################
# 
# this section (1) assesses input files to ensure they meet program input requirements and (2) provides a 
# user interaction allowing user to change input/output directories, input/output file names, and run the script
# with desired inputs.

# set workspace and ask if user wants to change it
theworkspace = "C:/overflights"
print("Overflights code 9-4-18. Workspace: "  )
print("")
print(theworkspace)
print("")
change = raw_input("Change workspace (y/n?)")

print("")

# if they do, ask for the new one.
if(change == "y"):
    theworkspace = raw_input("Enter new workspace: ")

    # if the user provided workspace is invalid, ask user again and again until valid path provided.
    while not os.path.exists(theworkspace):
        print("Invalid path.")
        print("")
        theworkspace = raw_input("Enter new workspace: ")
    # print the new workspace value
    print("Workspace changed to:")
    print(theworkspace)
    print("")
    
# ask if user wants to change the input point file.
change1 = raw_input("Input points: 2017trackpts.shp. Change?")
# default path:
pt_path = "2017trackpts.shp"
# if change is wanted, ask for new value and validate.
if( change1 == "y"):
    newpath = raw_input("New filename:")
        
    while not os.path.exists(theworkspace + "/" + newpath):
        newpath = raw_input("File does not exist. File must be a shapefile containing point data. Enter file name:")
        print("")
    # new path added here:
    pt_path = newpath

# ask if user wants to change the input lines file.
change2 = raw_input("Input lines: 2017tracklines.shp. Change?")
# default path
ln_path = "2017tracklines.shp"
# if change is wanted, ask for new value and validate.
if( change2 == "y"):
    newpath = raw_input("New filename:")
        
    while not os.path.exists(theworkspace + "/" + newpath):
        newpath = raw_input("File does not exist. File must be a shapefile containing line data. Enter file name:")
        print("")
    # new path added here
    ln_path = newpath

    
 # ask if user wants to change the input station file.   
change3 = raw_input("Input stations: 3d_stations.shp. Change?")
# default path
st_path = "3d_stations.shp"
# if change is wanted, ask for new value and validate.
if( change3 == "y"):
    newpath = raw_input("New filename:")
        
    while not os.path.exists(theworkspace + "/" + newpath):
        newpath = raw_input("File does not exist. File must be a shapefile containing station data. Enter file name:")
        print("")
    # new path added here:
    st_path = newpath    
    
    
    
# ask for user inputs: station code and search radius    
thecode = raw_input("Station code:")
thedistance = input("Search Radius (m):")
# station input file renamed below
thestations = st_path


# summarizes all input values prior to execution. 
print("")
print("Summary of analysis:")
print("")
print("Working directory: " + theworkspace)
print("Input point file: " + pt_path)
print("Input line file: " + ln_path)
print("Input sound station location file: 3d_stations.shp")
print("Search radius (m): " + str(thedistance))
tocontinue = raw_input("Press enter to begin event extraction.")
print("")

# set workspace and turn off overwrite protection
arcpy.env.workspace = theworkspace
arcpy.env.overwriteOutput = "True"

#####################################################
#################################################
####### Diagnose and report errors / problems with input files (*.shp)
###############################################

# make sure stations file has required field "code"
if len(arcpy.ListFields(thestations,"code"))==0:  
     print( "Missing field in station file: code  (see documentation)") 

# check submitted station code to ensure a matching value
station = arcpy.Select_analysis(thestations,"in_memory//station",'"code"= \'' + "UWBT" + '\'')

if(int(arcpy.GetCount_management(station).getOutput(0)) != 1):
    print("No matching stations. Please check input code. \n")

# check stations file to ensure it contains z data and point features
describe = arcpy.Describe(thestations)
if(describe.shapeType != u'Point'):
    print("Station file does not contain point features. (see documentation)")
if(describe.hasz != True):
    print("Station file does not contain z data. (see documentation)")
    
#####
# check input point file to ensure points are all present

# is it a point file
describe = arcpy.Describe(pt_path)
geometryType = describe.shapeType

if geometryType != u'Point':
    print("Point file does not contain point features. (see documentation)")

# does it contain the needed fields
if len(arcpy.ListFields(pt_path,"desc_"))==0:  
     print( "Missing field in point file: desc_  (see documentation)")  
if len(arcpy.ListFields(pt_path,"ltime"))==0:  
     print( "Missing field in point file: ltime  (see documentation)")  
if len(arcpy.ListFields(pt_path,"altitude"))==0:  
     print( "Missing field in point file: altitude  (see documentation)")      
        
#####        
# check input line file to ensure lines are in it and correct fields present

# is it a line file
describe = arcpy.Describe(ln_path)
geometryType = describe.shapeType
if geometryType != u'Polyline':
    print("Line file does not contain polyline features. (see documentation)")

# does it contain the needed fields
if len(arcpy.ListFields(ln_path,"desc_"))==0:  
     print( "Missing field in line file: desc_  (see documentation)")  
if len(arcpy.ListFields(ln_path,"ltime"))==0:  
     print( "Missing field in line file: ltime  (see documentation)")  
if len(arcpy.ListFields(ln_path,"altitude"))==0:  
     print( "Missing field in line file: altitude  (see documentation)") 


###############################################
###############################################


# run the code with the two provided inputs
extract_events(thecode,thedistance,flight_pts = pt_path,flight_lines = "2017tracklines.shp", outputdir = theworkspace,allstations = thestations )   

print("")
print("Output file saved to " + theworkspace + "/" + thecode + '_overflights.csv')
print("")



# In[2]:


# single_point = "moosecreekstation.shp"
# arcpy.ListFields(single_point)

# fields = arcpy.ListFields(single_point)

# for field in fields:
#     print("{0} is a type of {1} with a length of {2}"
#           .format(field.name, field.type, field.length))


# In[3]:


# # extract x coordinates of the two points
# xcoord = []
# with arcpy.da.SearchCursor(single_point, ['SHAPE@X']) as cursor:
#     for row in cursor:
#         xcoord.append(row[0]) 
        
# ycoord = []
# with arcpy.da.SearchCursor(single_point, ['SHAPE@Y']) as cursor:
#     for row in cursor:
#         ycoord.append(row[0]) 
        
# str(xcoord[0])

