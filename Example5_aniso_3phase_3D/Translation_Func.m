function f = Translation_Func(x)
x = normcdf(x);
f = x;
f(x<=0.3) = 0;
f(x>0.3 & x<=0.7) = 0.5;
f(x>0.7) = 1;


end