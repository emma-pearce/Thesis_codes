function [ro] = zdensity(v,z,Vice)


%This function calculates the density against depth from given velocities
%Version 1.0
%Julian Scott - British Antarctic Survey
%11 February 2009

%Density of ice
roi = 917; % in kgm^-3


%Density calculation from Kohnen (1972)
ro = roi .* (1+((Vice - v)./2.25).^1.22).^(-1);
ro = real(ro);


end