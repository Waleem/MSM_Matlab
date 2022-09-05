function g_m = g_test(m0,kbar)
    
    m1=2-m0;
    kbar2 = 2^kbar;
    g_m1 = [0:(kbar2-1)];

    for i = 1:(kbar2)
        g=1;
        for j = 0:(kbar-1)       
            if(bitand(g_m1(i),(2^j))~=0)    %
                g=g*m1;
            else g=g*m0;
            end
        end
        g_m1(i)=g;
    end
    
    g_m=sqrt(g_m1);
end