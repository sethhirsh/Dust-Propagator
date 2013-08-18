%    Copyright (c) 2010-2013, Delft University of Technology
%    All rights reserved.
%
%    Redistribution and use in source and binary forms, with or without modification, are
%    permitted provided that the following conditions are met:
%      - Redistributions of source code must retain the above copyright notice, this list of
%        conditions and the following disclaimer.
%      - Redistributions in binary form must reproduce the above copyright notice, this list of
%        conditions and the following disclaimer in the documentation and/or other materials
%        provided with the distribution.
%      - Neither the name of the Delft University of Technology nor the names of its contributors
%        may be used to endorse or promote products derived from this software without specific
%        prior written permission.
%
%    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
%    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
%    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
%    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
%    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
%    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
%    OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    Changelog
%      YYMMDD    Author            Comment
%
%
%    References
%
%    Notes
%


%This file plots the data for test cases between a 

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input deck.

% Set simulation data files.
dataFileDirectory = '/Users/sethmichaelhirsh/Desktop/Tudat/RawData/Cases/Correct Data/';

numCases = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loop through cases recording the data and the difference in the data in
%cell arrays.

for i = 1:numCases
    
    %Store all of COE data in cell array with numCases elements
    COEDataFiles{i} = csvread(strcat(dataFileDirectory,...
                                            'dustPropagationHistoryCOECase',... 
                                            num2str(i),...
                                            '_60s.dat'));
                                        
%Store all of the DNI data in cell array with numCases elements
    DNIDataFiles{i} = csvread(strcat(dataFileDirectory,...
                                            'dustPropagationHistoryDNICase',... 
                                            num2str(i),...
                                            '_60s.dat'));

%Record the difference in the DNI and COE data in a cell array with
%numCases elements
    DifferenceInDataFiles{i} = DNIDataFiles{i} - COEDataFiles{i};

end


%time data
%Note: Assumes that the time data is the same for all of the cases
time_data60s = COEDataFiles{1}(:,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following codes plots a set of subplots
%Each subplot is a plot of the difference of the DNI and COE integration
%methods for each orbital element and for each of the test cases with respect
%to time




%gets the current screen size and produces a figure with the specified
%position
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3)/2 scrsz(4)]);



%plot the Semi Major axis
semiMajorAxisPlot = subplot(3,2,1);
hold all
for i = 1:numCases
    plot(time_data60s, DifferenceInDataFiles{i}(:,2)./COEDataFiles{i}(:,2))
end
title('Difference in Semi-major Axis')
xlabel('Time (sec)')
ylabel('Semi-Major Axis (m)')
legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');



%Eccentricity
eccentricityPlot = subplot(3,2,2);
hold all
for i = 1:numCases
    plot(time_data60s, DifferenceInDataFiles{i}(:,3))
end
title('Difference in Eccentricity')
xlabel('Time (sec)')
ylabel('Eccentricity')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Inclination
inclinationPlot = subplot(3,2,3);
hold all
for i = 1:numCases
    plot(time_data60s, DifferenceInDataFiles{i}(:,4))
end
title('Difference in Inclination')
xlabel('Time (sec)')
ylabel('Inclination')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');


%Argument of Periapsis
argumentOfPeriapsisPlot = subplot(3,2,4);
hold all
for i = 1:numCases
    plot(time_data60s, DifferenceInDataFiles{i}(:,5))
end
title('Difference in Argument of Periapsis')
xlabel('Time (sec)')
ylabel('Argument of Periapsis')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Longitude of Ascending Node
longitudeOfAscendingNodePlot = subplot(3,2,5);
hold all
for i = 1:numCases
    plot(time_data60s, DifferenceInDataFiles{i}(:,6))
end
title('Difference in Longitude of Ascending Node')
xlabel('Time (sec)')
ylabel('Longitude of Ascending Node')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Longitude of Mean Anomaly
meanAnomalyPlot = subplot(3,2,6);
hold all
for i = 1:numCases
    
    %In the mean anomaly data there is an issue that in one data file the
    %mean anomaly may be recorded as very close to 360 degrees while in the
    %corresponding file the mean anomaly is recorded as very close to 0
    %degrees. Using the smallest magnitude function this code attempts to
    %handle the fact that the difference in these values would be 0 degrees
    %rather than -360 degrees
    meanAnomalyValues = smallestMagnitude(DifferenceInDataFiles{i}(:,7),360 ...
                        + DifferenceInDataFiles{i}(:,7));
                        
    plot(time_data60s,meanAnomalyValues)
end
title('Difference in Mean Anomaly')
xlabel('Time (sec)')
ylabel('Mean Anomaly')



% Set the 'super title' to appear above all of the subplots
suptitle('Difference in Orbital Elements with Step Size of 60 s')


%Legends are diplayed at the end because of interference from suptitle
legend(semiMajorAxisPlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
legend(eccentricityPlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
legend(inclinationPlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
legend(longitudeOfAscendingNodePlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
legend(argumentOfPeriapsisPlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
legend(meanAnomalyPlot,'Case I', 'Case II', 'Case III', 'Case IV', 'Case V');
