function [ue] = RichExtrap(u1, u2, u3)

k = 1 / ((u2-u3)/(u1-u2));

ue = u3 - (1/(k-1)) * (u2-u3);