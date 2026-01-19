% Homework 7, cant be done in python 

% simulation timing parameters
startTime = datetime("now");
stopTime = startTime + days(5);     % simulation length
sampleTime = 60; % seconds          % how often the simulation calculates data

sc = satelliteScenario(startTime, stopTime, sampleTime); % initialize the scenario with timing parameters

% orbital elements for the mission
rE = 6378;   % radius of earth [km]
z = 450;     % altitude above  [km]
ecc = 0.5;     % eccentricity
rp = rE + z;     % periapsis at desired altitude
a = rp/(1-ecc);  % semi major axis  [km]
i = 51.6;    % inclination     [deg]
Omega = 0;   % RAAN            [deg]
omega = 0;   % argument of periapsis  [deg]
theta = 0;   % true anomaly    [deg]

% create a single satellite orbit
v = satelliteScenarioViewer(sc, "PlaybackSpeedMultiplier", 50, "CameraReferenceFrame", "Inertial", "ShowDetails", false);

% setup the camera default position (optional really, it will be useful to move things)
latitude = 0;                     % degrees
longitude = 0;                    % degrees
height = rp*2*1000;                % meters
campos(v, latitude, longitude, height + 2000);

% create the satellite using the propagation engine
%sat = satellite(sc, a*1000 , ecc, i, Omega, omega, theta, "OrbitPropagator", "two-body-keplerian", "Viewer", v, "Name", "Sat310");
sat = satellite(sc,semiMajorAxis, ...
    eccentricity, ...
    inclination, ...
    rightAscensionOfAscendingNode, ...
    argumentOfPeriapsis, ...
    trueAnomaly, ...
    OrbitPropagator="numerical");
% Set the gravitational potential model used by the numerical propagator 
% to "oblate-ellipsode" (J2 perturbation).
numericalPropagator(sc,GravitationalPotentialModel="oblate-ellipsoid");
% draw the ground tracks for 48 hours forward (lead) and back (trail)
leadTime = 2*24*3600; % seconds
trailTime = leadTime*1000;
gt = groundTrack(sat, "LeadTime", leadTime, "TrailTime", trailTime);
o = orbit(sat, "LeadTime", leadTime, "TrailTime", trailTime);

% open the simulation
play(sc);