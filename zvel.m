function [V,Z] = zvel(a)


%This function calculates velocity against depth from the shallow
%refraction data
%Version 1.0
%Julian Scott - British Antarctic Survey
%11 February 2009

Distance = [1:30,30:10:1000,1000:100:10000]';
m = size(Distance,1);
for k = 1:m
    x(k) = Distance(k);
    recipvel(k) = a(1)*a(2)*exp(-a(2)*x(k)) + a(3)*a(4)*exp(-a(4)*x(k)) + a(5);
    velocity(k) = 1./recipvel(k);
    vel = velocity(k);
    z(k) = (1/pi) .* quad(@integrand, 0, x(k));
end

V = real(velocity);
Z = real(z);


function y = integrand(x)

%Version 1.0
%Julian Scott - British Antarctic Survey
%4 April 2008

s = a(1)*a(2)*exp(-a(2)*x) + a(3)*a(4)*exp(-a(4)*x) + a(5);
t = s .* vel;
y = acosh(t);


end

end