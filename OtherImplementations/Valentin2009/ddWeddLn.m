function  y = ddWeddLn(Ln_t, Ln_z, f, c_e)

if f == 1
    y = c_e * (1 + 3/(Ln_t^4*Ln_z^2));
elseif f == 2
    y = c_e *(1 + 3/(Ln_z^4*Ln_t^2));
elseif  f == 3
    y = c_e * 2 /(Ln_t*Ln_z)^3;
else
    exit('Wrong parameter for ddWddLn');
end