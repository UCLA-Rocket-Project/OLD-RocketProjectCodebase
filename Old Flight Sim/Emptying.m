function [d] = Emptying(initial, final, tbsteps, steps)
    d = linspace(initial, final, tbsteps);
    d(steps) = 0;