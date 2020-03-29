function [I] = PointMOI(componentMass,componentCG,rocketCG)
    d = rocketCG-componentCG;
    I = componentMass.*(d.^2);
end