s = 2:100; 
p_threshold = 1 - s.^(-1./s); 
[max_p, idx] = max(p_threshold); 

fprintf('The maximum prevalence p for which group testing is beneficial (T(s) < n) is p = %.8f\n', max_p);
