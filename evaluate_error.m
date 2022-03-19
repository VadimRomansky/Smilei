function [result] = evaluate_error(distribution, theta, beta, gamma, Np, startIndex, endIndex)
juttner(1:Np) = 0;
for i = 1:Np,
    juttner(i) = juttner_shifted_integrated(gamma(i),theta,beta);
end;
result = 0;
for i = startIndex:endIndex,
    result = result+((juttner(i) - distribution(i))^2)*(gamma(i)-gamma(i-1));
end;
end