function p_ = quick_val(x,alpha,beta,lambda)%  p_ = quick_val(x,alpha,beta,lambda) computes quick function at values%     in x given parameters alpha (threshold -- x at inflection point),%     beta (slope -- really shape -- term), and lambda (upper asymptote, %     or "lapse rate"). This version assumes a guessing rate of 0.5.% 	Returns vector or matrix, p, the same size as x. If x is a matrix, then % 		the function operates on each of the columns. alpha, beta & lambda must %		be	scalars.% see also quick_fit and quick_errp_ = 0.5 + (0.5 - lambda) * (1 - exp( -(x/alpha).^beta )); % quick%p_ = 0.5 + (0.5 - lambda) * (1 - 2.^( -(x/alpha).^beta )); % quick