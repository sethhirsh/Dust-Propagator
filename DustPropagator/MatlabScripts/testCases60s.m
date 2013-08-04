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

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input deck.

% Set simulation data files.
dataFileLocation = '/Users/sethmichaelhirsh/Desktop/Tudat/RawData/Cases/dustPropagationHistory'


% Case 1
dphCOE1_10s = csvread(strcat(dataFileLocation,'COECase1_10s.dat'));
dphCOE1_60s = csvread(strcat(dataFileLocation, 'COECase1_60s.dat'));
dphDNI1_10s = csvread(strcat(dataFileLocation, 'DNICase1_10s.dat'));
dphDNI1_60s = csvread(strcat(dataFileLocation, 'DNICase1_60s.dat'));

%Creates n x 6 arrays containing the difference in the orbital elements
diffOrbitalElem1_60s = dphDNI1_60s(:,2:end) - dphCOE1_60s(:,2:end);
diffOrbitalElem1_10s = dphDNI1_10s(:,2:end) - dphCOE1_10s(:,2:end);



% Case 2
dphCOE2_10s = csvread(strcat(dataFileLocation,'COECase2_10s.dat'));
dphCOE2_60s = csvread(strcat(dataFileLocation, 'COECase2_60s.dat'));
dphDNI2_10s = csvread(strcat(dataFileLocation, 'DNICase2_10s.dat'));
dphDNI2_60s = csvread(strcat(dataFileLocation, 'DNICase2_60s.dat'));

%Creates n x 6 arrays containing the difference in the orbital elements
diffOrbitalElem2_60s = dphDNI2_60s(:,2:end) - dphCOE2_60s(:,2:end);
diffOrbitalElem2_10s = dphDNI2_10s(:,2:end) - dphCOE2_10s(:,2:end);



% Case 3
dphCOE3_10s = csvread(strcat(dataFileLocation,'COECase3_10s.dat'));
dphCOE3_60s = csvread(strcat(dataFileLocation, 'COECase3_60s.dat'));
dphDNI3_10s = csvread(strcat(dataFileLocation, 'DNICase3_10s.dat'));
dphDNI3_60s = csvread(strcat(dataFileLocation, 'DNICase3_60s.dat'));

%Creates n x 6 arrays containing the difference in the orbital elements
diffOrbitalElem3_60s = dphDNI3_60s(:,2:end) - dphCOE3_60s(:,2:end);
diffOrbitalElem3_10s = dphDNI3_10s(:,2:end) - dphCOE3_10s(:,2:end);



% Case 4
dphCOE4_10s = csvread(strcat(dataFileLocation,'COECase4_10s.dat'));
dphCOE4_60s = csvread(strcat(dataFileLocation, 'COECase4_60s.dat'));
dphDNI4_10s = csvread(strcat(dataFileLocation, 'DNICase4_10s.dat'));
dphDNI4_60s = csvread(strcat(dataFileLocation, 'DNICase4_60s.dat'));

%Creates n x 6 arrays containing the difference in the orbital elements
diffOrbitalElem4_60s = dphDNI4_60s(:,2:end) - dphCOE4_60s(:,2:end);
diffOrbitalElem4_10s = dphDNI4_10s(:,2:end) - dphCOE4_10s(:,2:end);



% Case 5
dphCOE5_10s = csvread(strcat(dataFileLocation,'COECase5_10s.dat'));
dphCOE5_60s = csvread(strcat(dataFileLocation, 'COECase5_60s.dat'));
dphDNI5_10s = csvread(strcat(dataFileLocation, 'DNICase5_10s.dat'));
dphDNI5_60s = csvread(strcat(dataFileLocation, 'DNICase5_60s.dat'));

%Creates n x 6 arrays containing the difference in the orbital elements
diffOrbitalElem5_60s = dphDNI5_60s(:,2:end) - dphCOE5_60s(:,2:end);
diffOrbitalElem5_10s = dphDNI5_10s(:,2:end) - dphCOE5_10s(:,2:end);



%time data
time_data60s = dphCOE1_60s(:,1);
time_data10s = dphCOE1_10s(:,1);



fig1 = figure;


%Semi Major axis
subplot(3,2,1)
%axis([0.0 7e5 -200 50])
hold all
plot(time_data60s,diffOrbitalElem1_60s(:,1))
plot(time_data60s,diffOrbitalElem2_60s(:,1))
plot(time_data60s,diffOrbitalElem3_60s(:,1))
plot(time_data60s,diffOrbitalElem4_60s(:,1))
plot(time_data60s,diffOrbitalElem5_60s(:,1))

title('Difference in Semi-major Axis')
xlabel('Time (sec)')
ylabel('Semi-Major Axis')


legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Eccentricity
subplot(3,2,2)
hold all
plot(time_data60s,diffOrbitalElem1_60s(:,2))
plot(time_data60s,diffOrbitalElem2_60s(:,2))
plot(time_data60s,diffOrbitalElem3_60s(:,2))
plot(time_data60s,diffOrbitalElem4_60s(:,2))
plot(time_data60s,diffOrbitalElem5_60s(:,2))

title('Difference in Eccentricity')
xlabel('Time (sec)')
ylabel('Eccentricity')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Inclination
subplot(3,2,3)
hold all
plot(time_data60s,diffOrbitalElem1_60s(:,3))
plot(time_data60s,diffOrbitalElem2_60s(:,3))
plot(time_data60s,diffOrbitalElem3_60s(:,3))
plot(time_data60s,diffOrbitalElem4_60s(:,3))
plot(time_data60s,diffOrbitalElem5_60s(:,3))

title('Difference in Eccentricity')
xlabel('Time (sec)')
ylabel('Eccentricity')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');


%Argument of Periapsis
subplot(3,2,4)
hold all
plot(time_data60s,diffOrbitalElem1_60s(:,4))
plot(time_data60s,diffOrbitalElem2_60s(:,4))
plot(time_data60s,diffOrbitalElem3_60s(:,4))
plot(time_data60s,diffOrbitalElem4_60s(:,4))
plot(time_data60s,diffOrbitalElem5_60s(:,4))

title('Difference in Argument of Periapsis')
xlabel('Time (sec)')
ylabel('Argument of Periapsis')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Longitude of Ascending Node
subplot(3,2,5)
hold all
plot(time_data60s,diffOrbitalElem1_60s(:,5))
plot(time_data60s,diffOrbitalElem2_60s(:,5))
plot(time_data60s,diffOrbitalElem3_60s(:,5))
plot(time_data60s,diffOrbitalElem4_60s(:,5))
plot(time_data60s,diffOrbitalElem5_60s(:,5))

title('Difference in Longitude of Ascending Node')
xlabel('Time (sec)')
ylabel('Longitude of Ascending Node')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

%Longitude of Mean Anomaly
subplot(3,2,6)
hold all
diffMeanAnomaly1 = smallestMagnitude(diffOrbitalElem1_60s(:,6),(360 + diffOrbitalElem1_60s(:,6)));
diffMeanAnomaly2 = smallestMagnitude(diffOrbitalElem2_60s(:,6),(360 + diffOrbitalElem2_60s(:,6)));
diffMeanAnomaly3 = smallestMagnitude(diffOrbitalElem3_60s(:,6),(360 + diffOrbitalElem3_60s(:,6)));
diffMeanAnomaly4 = smallestMagnitude(diffOrbitalElem4_60s(:,6),(360 + diffOrbitalElem4_60s(:,6)));
diffMeanAnomaly5 = smallestMagnitude(diffOrbitalElem5_60s(:,6),(360 + diffOrbitalElem5_60s(:,6)));


plot(time_data60s,diffMeanAnomaly1);
plot(time_data60s,diffMeanAnomaly2);
plot(time_data60s,diffMeanAnomaly3);
plot(time_data60s,diffMeanAnomaly4);
plot(time_data60s,diffMeanAnomaly5);

title('Difference in Mean Anomaly')
xlabel('Time (sec)')
ylabel('Mean Anomaly')

legend('Case I', 'Case II', 'Case III', 'Case IV', 'Case V');

suptitle('Difference in Orbital Elements with Step Size of 60 s')















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read and store simulation data files.
% First column is epoch, subsequent columns are Cartesian state elements.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot orbits of Asterix and Obelix.
% Simulation data is given in m and m/s, so needs to be converted to km.



