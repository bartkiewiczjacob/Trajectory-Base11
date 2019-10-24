function [g] = gravity(alt)

G = 6.673E-11;
Mearth = 5.98E24;

g = G*Mearth / alt^2;

