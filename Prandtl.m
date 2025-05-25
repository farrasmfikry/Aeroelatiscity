function F = Prandtl(B,R,r,phi)
    F = (2/pi)*acos(exp((-B/2)*((R-r)/(r*sin(abs(phi))))));
end
