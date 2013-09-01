#######################################################################

# Setup input deck

#########################################################################

# Set input directory
inputDirectory = ""

# Set case name
caseName = "KeplerOrbitTest"
 
# Set number of Cases
numCases = 5

# Set stepsize [seconds]
stepsize = 60

#Set absolute path directory
outputPath = ""

#Set figure DPI
figureDPI = 800

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


#Create dictionary to store all inputed data from the COE directories
COEDataDictionary = {}

#Create dictionary to store all inputed data from the DNI directories
DNIDataDictionary = {}


#Loop through the cases and store each case in its corresponding dictionary
for i in range(numCases):

	caseNumber = str(i + 1)

	# Read COE data file and convert it to a list
	# A list containing a sublist with 7 elements is expected
	COENameOfFile = 'dustStateHistoryKeplerElementsBenchmark_case' + caseNumber + '_' + str(stepsize) + 's.dat'
	COEDataList = list(reader(open(inputDirectory + COENameOfFile,'rt')))


	# Create an array of the orbital elements (excluding time which should be in the first column)
	COEOrbitalElementArray = np.array(COEDataList, dtype = float)[:,1:]

	# Store time intervals in separate array
	COETimeData = np.array(COEDataList)[:,0]

	# Store array of orbital elements in in dictionary
	COEDataDictionary['Case' + caseNumber + str(stepsize) + 's'] = COEOrbitalElementArray

	# Read DNI data file and convert it to a list
	# A list containing a sublist with 7 elements is expected
	DNINameOfFile = 'dustStateHistoryKeplerElements_case' + caseNumber + '_' + str(stepsize) + 's.dat'
	DNIDataList = list(reader(open(inputDirectory + DNINameOfFile,'rt')))

	# Create an array of the orbital elements from the DNI data list
	# The first column which stores the epochs is excluded
	DNIOrbitalElementArray = np.array(DNIDataList, dtype = float)[:,1:]

	# Store time intervals in separate array
	DNITimeData = np.array(DNIDataList)[:,0]

	# Story array of orbital elements in in dictionary
	DNIDataDictionary['Case' + caseNumber + str(stepsize) + 's'] = DNIOrbitalElementArray

	# Check to make sure that COE and DNI time data is the same
	assert all(COETimeData == DNITimeData), 'COE and DNI do not have the same time data'


###########################################################################################

# Define functions that can be used to calculate error
# Other functions to manipulate data can be defined here

########################################################################################################

def absoluteError(testArray, benchmarkArray):
	# Recommended input parameters are numpy arrays
	assert type(testArray) == np.ndarray and type(benchmarkArray) == np.ndarray, 'Input should be numpy arrays'

	# Return absolute error of test array with respect to benchmark array
	return testArray - benchmarkArray

def relativeError(testArray, benchmarkArray):
	# Recommended input parameters are numpy arrays
	assert type(testArray) == np.ndarray and type(benchmarkArray) == np.ndarray, 'Input should be numpy arrays'

	# Return relative error of test array with respect to benchmark array
	return (testArray - benchmarkArray) / benchmarkArray


####################################################################################################

# Create a list where each element is a list containing:
# The name of an orbital element
# The unit of measurement for that orbital element
# The function that will be used to calculate the error in the orbital element data
# The name of the function described on the previous line

orbitalElements = [ \
					['Semi-major Axis', 'meters', absoluteError, 'Absolute Difference'], \
					['Eccentricity', None, absoluteError, 'Absolute Difference'], \
					['Inclination' , 'degrees' , absoluteError, 'Absolute Difference'], \
					['Argument of Periapsis' , 'degrees', absoluteError, 'Absolute Difference'], \
					['Longitude of Ascending Node' , 'degrees', absoluteError, 'Absolute Difference'], \
					['Mean Anomaly' , 'degrees', absoluteError, 'Absolute Difference'] \
					]

####################################################################################################

# Plot the difference between the DNI and COE integrations

######################################################################################################

# Define a function to return an array of the values of smallest magnitude from two arrays
# This function will be used to edit the mean anomaly data
def smallestMagnitude(npArray1, npArray2):
	# Expected input are two one-dimensional numpy arrays of equal size

	# Check that both numpy arrays are of equal length
	assert len(npArray1) == len(npArray2) , 'Numpy arrays are not equal sizes'

	# Create new numpy array to contain the value with the largest magnitude for each
	# set of corresponding elements of the array
	smallestMagnitudeArray = np.zeros(len(npArray1), dtype = float)

	# Loop through the two arrays and put the value of smallest magnitude in an array 
	for i in range(len(npArray1)):

		# Check which element is smaller in magnitude and store it in smallestMagnitudeArray
		if abs(npArray1[i]) < abs(npArray2[i]):
			smallestMagnitudeArray[i] = npArray1[i]
		else:
			smallestMagnitudeArray[i] = npArray2[i]

	# Return the array of values with the smallest magnitude
	return smallestMagnitudeArray

#######################################################################################################


# Create a figure
plt.figure()

# Add a super title
#plt.suptitle('Absolute Difference between DNI and COE Integration')

# Loop through the orbital elements and plot their differences for each case
for i in range(len(orbitalElements)):
	
	# Create subplot for each orbital element
	ax = plt.subplot(3,2,i+1)

	#Get the orbital element name
	orbitalElementName = orbitalElements[i][0]

	#Get the unit of the orbital element if it exists
	orbitalElementUnit = orbitalElements[i][1]

	# Get the function that will be used to calculate the error for the orbital element
	errorFunction = orbitalElements[i][2]

	# Get the name of the type of error calculation
	errorType = orbitalElements[i][3]

	# Create list to hold all of the plots (will be used to create legend)
	casePlotList = []

	# Create list to hold labels (for legend) for each case
	legendLabels = []

	# Loop through each case and plot the error in the data
	for j in range(numCases):

		# Case numbers start at 1
		caseNumber = str(j + 1)

		# String which is key for dictionary
		currentKey = 'Case' + caseNumber + str(stepsize) + 's'

		# Calculate error in DNI data with respect to COE data
		errorData = errorFunction(DNIDataDictionary[currentKey][:,i],COEDataDictionary[currentKey][:,i])

		# Alter error data when using the mean anomaly
		# Issue stems from the fact that 0 degrees and 360 degrees are equivalent
		if orbitalElementName == 'Mean Anomaly':

			errorData = smallestMagnitude(errorData, errorData + 360.0 )

		# Plot the absolute difference between the DNI and COE data vs. time
		casePlot, = ax.plot(COETimeData, errorData, label = caseNumber)

		# Append the casePlot to a list (this is necessary to create the legend)
		casePlotList.append(casePlot)

		# Append label to list (for legend)
		legendLabels.append('Case ' + caseNumber)

	#Add title and labels to subplot

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
	plt.title(errorType + ' in ' + orbitalElementName, fontsize = 8)

	# Create legend
	plt.legend(casePlotList,legendLabels, fontsize = 6)

	# Set linewidth
	pylab.setp(casePlotList, linewidth = 0.5)

# Prevent overlapping of subplots and title
plt.tight_layout(h_pad = 2.0)

# Prevents Overlapping of supertitle
#plt.subplots_adjust(top = 0.90)

# Save file 
plt.savefig(outputPath + 'DNIErrorInOrbitalElements' + str(stepsize) + 's.png' , \
            dpi = figureDPI)   


# Close plot 
plt.close()