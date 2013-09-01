#######################################################################

# Setup input deck

#########################################################################

# Set input directory
inputDataFile = ""

# Set case name
caseName = ""
 
# Set stepsize [seconds]
stepsize = 10

#Set absolute path directory
outputPath = ""

#Set figure DPI
figureDPI = 800

#Set the super title for the plot
supertitle = ""

###########################################################################################

#Import Python packages and modules

#################################################################################################

import matplotlib.pyplot as plt
import os
from csv import reader
import numpy as np
import pylab


###################################################################################################

# Load input data and record it in dictionaries

###################################################################################################

#Check if output directory exists; if not, create it
if not os.path.exists(outputPath):
	os.makedirs(outputPath)


# Convert Data File To a list
dataList = list(reader(open(inputDirectory + inputDataFile,'rt')))


# Create an array of the orbital elements (excluding time which should be in the first column)
orbitalElementArray = np.array(COEDataList, dtype = float)[:,1:]

# Store time intervals in separate array
timeData = np.array(DataList)[:,0]


###########################################################################################

# Create a list where each element is a list containing:
# The name of an orbital element
# The unit of measurement for that orbital element

###########################################################################################

orbitalElements = [ \
					['Semi-major Axis', 'meters'], \
					['Eccentricity', None], \
					['Inclination' , 'degrees'], \
					['Argument of Periapsis' , 'degrees'], \
					['Longitude of Ascending Node'], \
					['Mean Anomaly' , 'degrees'] \
					]

####################################################################################################

# Plot the data

######################################################################################################


# Create a figure
plt.figure()

# Add a super title
plt.suptitle(supertitle)

# Loop through the orbital elements and plot their differences for each case
for i in range(len(orbitalElements)):
	
	# Create subplot for each orbital element
	ax = plt.subplot(3,2,i+1)

	# Get the orbital element name
	orbitalElementName = orbitalElements[i][0]

	# Get the unit of the orbital element if it exists
	orbitalElementUnit = orbitalElements[i][1]

	# Get the data for the orbital element
	orbitalElementData = OrbitalElementArray[:,i]

	# Plot the data
	dataPlot = plt.plot(timeData, orbitalElementData)

	# Set linewith
	pylab.setp(dataPlot, linewidth = 10)

	# Sets default font size
	plt.rc('font', **{'size':'6'})

	# Format tick marks on x-axis
	plt.xticks(fontsize = 6)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(-3,3))

	# Format tick marks on y-axis
	plt.yticks(fontsize = 6)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))

	# Add label to x-axis
	plt.xlabel('Time [seconds]', fontsize = 7)

	# Add label to y-axis
	# Adds the orbital elements units to the label if it has a unit
	if orbitalElementUnit is None or errorFunction == relativeError:
		plt.ylabel(orbitalElementName, fontsize = 7)
	else:
		plt.ylabel(orbitalElementName + '\n' + '[' + orbitalElementUnit + ']', fontsize = 8)

	plt.title(orbitalElementName, fontsize = 8)

# Prevent overlapping of subplots and title
plt.tight_layout(h_pad = 2.0)

# Prevents Overlapping of supertitle
plt.subplots_adjust(top = 0.90)

# Save file 
plt.savefig(outputPath + caseName + str(stepsize) + 's.png' , \
            dpi = figureDPI)   


# Close plot 
plt.close()