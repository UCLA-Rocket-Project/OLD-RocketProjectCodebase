function [I] = CylinderMOI(mass, outradius, inradius, length)
I = (mass/12).*(3*(outradius^2 + inradius^2) + length.^2);