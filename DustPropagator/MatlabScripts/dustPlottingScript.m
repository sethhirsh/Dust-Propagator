%    Copyright (c) 2010-2013, S. Hirsh
%    Copyright (c) 2010-2013, Delft University of Technology
%    All rights reserved.
%    See COPYING for license details.
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

% Set Asterix simulation data file
dustSimulationDataFile = '/Users/sethmichaelhirsh/Desktop/Tudat/dustPropagationHistory.dat'
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Sun.
% This code was provided by J. Melman, TU Delft.
% Axes are given in km.

% Load data file.
load( 'topo.mat' )

% Set radius of Sun in km.
%Note: the Sun appears to be very small in respect to the size of 1 AU
r_S = 695500;


% Create geometric sphere.
[x y z] = sphere( 50 );

% Set plot properties.
hold on
colormap( hot );

props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Sun.
surface( r_S * x, r_S * y, r_S * z, props );
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read and store simulation data files.
% First column is epoch, subsequent columns are Cartesian state elements.

% Read and store dust simulation data
dustSimulationData = csvread( dustSimulationDataFile );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot orbits of dust particle
% Simulation data is given in m and m/s, so needs to be converted to km.

% Convert data to km.
dustSimulationData = dustSimulationData / 1000;

% Plot orbits.
grid on;
xlabel( 'Cartesian x-position [km]' );
ylabel( 'Cartesian y-position [km]' );
zlabel( 'Cartesian z-position [km]' );

plot3( dustSimulationData(:,2),dustSimulationData(:,3),...
    dustSimulationData(:,4), 'LineWidth', 3 );


for i = 1:size(dustSimulationData,1)
    
    % Plot instantaneous position of dust particle
    figureHandleDust = plot3( dustSimulationData(i,2),...
                                 dustSimulationData(i,3),...
                                 dustSimulationData(i,4), 'o',...
                                 'MarkerFaceColor', 'y',...
                                 'MarkerEdgeColor', 'k',...
                                 'MarkerSize', 10 );
 
    
   % Capture frame for movie.
    simulationMovieFrames(i) = getframe;
    
    delete( figureHandleDust );
  
end

% Play movie.
movie(simulationMovieFrames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
