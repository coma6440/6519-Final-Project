function [a, e, i, OMEGA, omega, theta, h] = get_orbital_elements(R, V, mu)
% get_orbital_elements
%
%   Function to calculate the orbital elements of a satellite based on the
%   position and velocity vectors of the spacecraft at a given time t1
%
%   Inputs:
%       R - [3x1] [km] Position vector in inertial coordinates of the
%           spacecraft at time t1
%       V - [3x1] [km/s] Velocity vector in inertial coordinates of the
%           spacecraft at time t1
%       mu - [km^3/s^2] Gravitational Parameter of the Central Body
%
%   Outputs:
%       a - [km] Semi-Major Axis
%       e - Eccentricity
%       i - [deg] Inclination
%       OMEGA - [deg] Right Ascension of the Ascending Node
%       omega - [deg] Argument of Periapsis
%       theta - [deg] True Anomaly

r = norm(R); % [km] Distance from central body
v = norm(V); % [km/s] Speed

h = cross(R, V); % [km^2/s] Specific Angular Momentum
E = v^2/2 - mu/r; % [km^2/s^2] Specific Energy

a = -mu/(2*E); % [km] Semi-Major Axis
e = sqrt(1 + 2*norm(h)^2*E/mu^2); % Eccentricity
i = acosd(h(3)/norm(h)); % [deg] Inclination

n = cross([0; 0; 1], h); % Line of Nodes
OMEGA = wrapTo180(acosd(dot(n, [1; 0; 0])/norm(n))); % [deg] Right Ascention of the Ascending Node

if dot(n, [0; 1; 0]) > 0
    OMEGA = abs(OMEGA);
else
    OMEGA = -abs(OMEGA);
end

e_vec = cross(V, h)/mu - R/r; % Eccentricity Vector

omega = wrapTo180(acosd(dot(n, e_vec)/(norm(n)*e))); % Argument of Periapsis

if dot(e_vec, [0; 0; 1]) > 0
    omega = abs(omega);
else
    omega = -abs(omega);
end

theta = wrapTo180(acosd(dot(R, e_vec)/(r*e))); % True Anomaly

if dot(R, V) > 0
    theta = abs(theta);
else
    theta = -abs(theta);
end


end