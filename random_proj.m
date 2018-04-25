% measure the expected value of a random projection
% analysis for the paper by Dasgupta1999

% uniform projection of a circle to 1 d
theta = rand(100000,1)*2*pi;

L = abs(cos(theta));

mean(L)
display(2/pi);

% restrict thinking to only 1 quadrant. all other quandrants through
% symmetry
% for d = 2, k = 1, sample \theta uniformly and E[L] = int_[0,pi/2] cos t
% dt/(pi/2) = 2/pi
% for d = 3, k = 2, p(theta,phi) = [1/(pi/2), sin \phi].
% E[L] = pi/4 


% uniform projection of d to k 
clear all;
d = 3;
k = 2;

X = randn(10000,d);
normX = sqrt(sum(X.^2,2) );
X_d = diag(normX)\X ; 
clear X;
%L_d = sqrt( sum(X_d.^2,2) );
%hist(L_d);
X_k = X_d(:,1:k); 
clear X_d
L_k = sqrt( sum(X_k.^2,2) );
hist(L_k);
display(mean(L_k));
display(pi/4);


%Anyhoo - i don't see how dasgupta gets EL = k/d