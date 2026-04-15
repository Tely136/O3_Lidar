function fl = gen_framelength(M1,M2,h1,h2,h)
    N1 = (M1-1)/2;
    N2 = (M2-1)/2;
    
    fl = 2.*(round(((N2-N1)/(h2-h1)).*(h-h1) + N1))+1;
    fl(h < h1) = M1;
    fl(h > h2) = M2;
end