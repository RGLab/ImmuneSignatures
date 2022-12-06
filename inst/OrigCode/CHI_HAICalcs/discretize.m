function r=discretize(x, top, bot )
    assert(top>0 & top <100 );
    assert(100-top > bot );
    assert (bot>0 & bot<100);
    
    top_cut = prctile(x, 100-top);
    low_cut = prctile(x, bot);
    
    r = ones(size(x));
    r( x<=low_cut ) = 0;
    r( x>=top_cut ) = 2;
    
    r( isnan(x)) = NaN;
end