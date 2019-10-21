function sol = hist_sizing(d,n)

df = d/12;
totA = (12.1*df -2.61);
finA = totA/n*144*.7; 
% This 0.7 factor corrects for the improvements in material it is derived
% based on trial and error.

a = 1.6*d;
b = 1.2*d;
h = 2*finA/2.8/d;

sol = [a,b,h];

end
