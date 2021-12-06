function [ztwtt,twtt] = twoway(v,z,Vice)


%This function calculates the two-way travel time against depth from given velocities
%Version 1.0
%Julian Scott - British Antarctic Survey
%11 February 2009

%First interpolate to give a velocity at depths of 0.5,1.5,2.5,... metres
%up to the deepest we have a velocity value for
[z,rows,columns] = unique(z); %need to make sure values are unique for interp
v = v(rows);

zmax = floor(z(end));
zi = 0.5:1:zmax-0.5;
vi = interp1(z,v,zi,'linear','extrap');

%Now calculate the travel time for each metre travelled from the surface using the
%velocity at the centre of each depth range

twtti = 2./vi; %For each step

%Check size and if less than 150 m extrapolate this down to 150 m using Vice
[n] = size(twtti,2);
if n < 150
    twtti(1,n+1:150) = 2./Vice;

    %Now calculate total cumulative travel time
    twtt(1) = twtti(1);
    for k = 2:150
    twtt(k) = twtti(k) + twtt(k-1); %Total culmulative
    end
    ztwtt = [1:150];
else
    %Now calculate total cumulative travel time
    twtt(1) = twtti(1);
    for k = 2:zmax
    twtt(k) = twtti(k) + twtt(k-1); %Total culmulative
    end
    ztwtt = [1:zmax];
end


end