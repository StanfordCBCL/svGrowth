function  y = dWmdLn(Ln, c_m2, c_m3, c_m2_c, c_m3_c)

if Ln >= 1
    y = c_m2*Ln*(Ln^2 - 1)*exp(c_m3*(Ln^2 - 1)^2);
else
    y= c_m2_c*Ln*(Ln^2-1.0)*exp(c_m3_c*(Ln^2-1.0)^2);
end