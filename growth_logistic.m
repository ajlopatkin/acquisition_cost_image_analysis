function output = growth_logistic(P,t)

output = P(1) ./ (1 + exp(4 * P(2) / P(1) .* (P(3) - t) + 2) );

end