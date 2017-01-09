function d = KLNormMulti(m1,S1,m2,S2)

if((sum(isnan(S1(:)))>0) || (sum(isnan(S2(:)))>0))
    d =NaN;
else
    d = 0.5*( trace(inv(S2)*S1) + (m2-m1)*inv(S2)*(m2-m1)' - length(m1) -log(det(S1)/det(S2)) );
end