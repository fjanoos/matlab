function [out1] = lor_eq(x0,y0,z0,beta,sigma,ro)
out1 = [sigma*(y0-x0);x0*(ro-z0)-y0;x0*y0-beta*z0];