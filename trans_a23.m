function [a_23] = trans_a23(theta_blade)

a_23 = [cos(theta_blade) sin(theta_blade) 0;-sin(theta_blade) cos(theta_blade) 0;0 0 1];

end