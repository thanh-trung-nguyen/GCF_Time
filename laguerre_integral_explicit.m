function f = laguerre_integral_explicit(n,x,x0)

% calculate the integral of the  Laguerre's function using explicit formulae given by Maple.
% Maximum dimension is 10.

if nargin < 3
    x0 = 0;
end
x = x - x0;

switch n
    case 0
        f = 2;
    case 1
        f = -2*(1+x);
    case 2
        f = 2 + x.^2;
    case 3
        f = -(6 + (6 + (-3 + x).*x).*x)/3;
    case 4
        f = (24 + (24 + (-8 + x).*x).*x.*x)/12;
    case 5
        f = -(120 + (120 + (-120 + (80 + (-15 + x).*x).*x).*x).*x)/60;
    case 6
        f = (720 + (1080 + (-720 + (210 + (-24 + x).*x).*x).*x).*x.*x)/360;
    case 7
        f = -(5040 + (5040 + (-7560 + (7560 + (-2730 + (462 + (-35 + x).*x).*x).*x).*x).*x).*x)/2520;
    case 8
        f = (40320 + (80640 + (-80640 + (36960 +(-8064 + (896 +(-48 + x).*x).*x).*x).*x).*x).*x.*x)/20160;
    case 9
        f = -(362880 + (362880 + (-725760 + (967680 +(-514080 + (139104 + (-20160 + (1584 +(-63 + x).*x).*x).*x).*x).*x).*x).*x).*x)/181440;
    case 10
        f = (1/1814400)*(3628800+9072000*x.^2-12096000*x.^3+7560000*x.^4-2419200*x.^5+433440*x.^6-44640*x.^7+2610*x.^8-80*x.^9+x.^10);
    case 11
        f = -(1/19958400)*(39916800+39916800*x-99792000*x.^2+166320000*x.^3-116424000*x.^4+43243200*x.^5-9203040*x.^6+1172160*x.^7-90090*x.^8+4070*x.^9-99*x.^10+x.^11);
    case 12
        f = (1/239500800)*(479001600+1437004800*x.^2-2395008000*x.^3+1896048000*x.^4-798336000*x.^5+196922880*x.^6-29842560*x.^7+2839320*x.^8-168960*x.^9+6072*x.^10-120*x.^11+x.^12);
    case 13
        f =-(1/3113510400)*(6227020800+6227020800*x-18681062400*x.^2+37362124800*x.^3-32432400000*x.^4+15308092800*x.^5-4289725440*x.^6+753667200*x.^7-85405320*x.^8+6297720*x.^9-298584*x.^10+8736*x.^11-143*x.^12+x.^13);
    case 14
        f =(1/43589145600)*(87178291200+305124019200*x.^2-610248038400*x.^3+584821036800*x.^4-305124019200*x.^5+95775039360*x.^6-19130791680*x.^7+2514592080*x.^8-221020800*x.^9+12996984*x.^10-502320*x.^11+12194*x.^12-168*x.^13+x.^14);
        
    case 15
        f =-(1/653837184000)*(1307674368000+1307674368000*x-4576860288000*x.^2+10679340672000*x.^3-11060745696000*x.^4+6331323398400*x.^5-2199435638400*x.^6+492194102400*x.^7-73589115600*x.^8+7506298800*x.^9-526485960*x.^10+25257960*x.^11-810810*x.^12+16590*x.^13-195*x.^14+x.^15);
    otherwise
        error('n must be from 0 to 15');
end
f = f.*exp(-x/2);
f = f*(-1)^n; %to make all the coefficient positives