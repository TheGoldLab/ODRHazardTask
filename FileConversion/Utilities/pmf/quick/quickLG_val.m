function p_ = quickLG_val(x,alpha,beta,lambda,gamma)%  p_ = quickLG_val(x,alpha,beta,lambda) computes quick function at values%     in x given parameters alpha (threshold -- x at inflection point),%     beta (slope -- really shape -- term), lambda (upper asymptote, %     or "lapse rate"), and gamma (lower asymptote, or "guess rate"). % 	Returns vector or matrix, p, the same size as x. If x is a matrix, then % 		the function operates on each of the columns. alpha, beta & lambda must %		be	scalars.% see also quickLG_fit and quickLG_errp_ = gamma + (gamma - lambda) * (1 - exp( -(x/alpha).^beta )); % quick