function [t, p, nd] = standard_atmosphere(z_m)
    t = NaN(size(z_m));
    p = NaN(size(z_m));
    
    g = 9.80665;
    M = 0.0289644;
    R = 8.3144598;
    NA = 6.02214076 * 10^23;
    
    for i = 1:length(z_m)
        zi = z_m(i);
    
        if zi >=0 && zi < 11000
            T0 = 288.15;
            p0 = 101325;
            z0 = 0;
            L = -0.0065;
    
        elseif zi >= 11000 && zi < 20000
            T0 = 216.65;
            p0 = 22632.1;
            z0 = 11000;
            L = 0;
    
    
        elseif zi >= 20000 && zi < 32000
            T0 = 216.65;
            p0 = 5474.89;
            z0 = 20000;
            L = 0.001;
    
        elseif zi >= 32000 && zi < 47000
            T0 = 228.65;
            p0 = 868.019;
            z0 = 32000;
            L = 0.0028;
    
        elseif zi >= 47000 && zi < 51000
            T0 = 270.65;
            110.906;
            z0 = 47000;
            L = 0;
    
        elseif zi >= 51000 && zi < 71000
            T0 = 270.65;
            p0 = 66.9389;
            z0 = 51000;
            L = -0.0028;
    
        elseif zi >= 71000 && zi < 84852
            T0 = 214.65;
            p0 = 3.95642;
            z0 = 71000;
            L = -0.002;
    
        elseif zi > 84852
            T0 = 186.946;
            p0 = 0;
            z0 = 84852;
            L = 0;
    
        else
            %  error
        end
        
        t_temp = T0 + L * (zi - z0);
    
        if L == 0
            p(i) = p0 * exp(-g*M*(zi-z0)/(R*T0));
        else
            p(i) = p0 * (T0/t_temp).^(g*M/R/L);
        end
    
        t(i) = t_temp;
    end

    nd = p .* NA ./ (R .* t); % m^-3
end