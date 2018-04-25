function [ ] = lorenz_pdf ( )
%This function will be used to see how badly the original pdf is skewed,
%this is a test done to see the effect of propagating a pdf through a
%nonlinear function. 

x = 0.01.*randn(1,10000) + 1;
y = 0.01.*randn(1,10000) + 0;
z = 0.01.*randn(1,10000) + 0;
beta  = 8/3; sigma = 10; ro = 28; h = 0.05; %Parameters for the Lorenz equation

%Now the initial points will be used as an input to Lorenz_att.m and a
%future point will be obtained. The duration of integration will be
%progressively increased so as to obtain enough data points.

%Defining the maximum iteration for the Lorenz_att.m and also initializing
%variables to construct the pdf.
max_iteration = 2000;
x_out = zeros(size(x,2),max_iteration);
y_out = zeros(size(x,2),max_iteration);
z_out = zeros(size(x,2),max_iteration);

[nx xout] = hist(x,100);
[ny yout] = hist(y,100);
[nz zout] = hist(z,100);

subplot(3,1,1)
    plot(xout,nx./size(x,2),'black');
    xlabel('X values from the Lorenz attractor')
    ylabel('Probability values')
    
    subplot(3,1,2)
    plot(yout,ny./size(y,2),'blue');
    xlabel('Y values from the Lorenz attractor')
    ylabel('Probability values')
    
    subplot(3,1,3)
    plot(zout,nz./size(z,2),'black');
    xlabel('Z values from the Lorenz attractor')
    ylabel('Probability values')
    
    cd('../images/');
    print -djpeg pdf_1
    close all
    cd('../monte_carlo/')

    matlabpool(3)
for iteration = 1:max_iteration
    parfor count = 1:size(x,2)
        [a,b,c] = lorenz_att(x(count),y(count),z(count),beta,sigma,ro,h,iteration);
        x_out(count,iteration) = a(end);
        y_out(count,iteration) = b(end);
        z_out(count,iteration) = c(end);
    end
    [nx xout] = hist(x_out(:,iteration),100);
    [ny yout] = hist(y_out(:,iteration),100);
    [nz zout] = hist(z_out(:,iteration),100);
    
    subplot(3,1,1)
    plot(xout,nx./size(x,2),'black');
    xlabel('X values from the Lorenz attractor')
    ylabel('Probability values')
    
    subplot(3,1,2)
    plot(yout,ny./size(y,2),'blue');
    xlabel('Y values from the Lorenz attractor')
    ylabel('Probability values')
    
    subplot(3,1,3)
    plot(zout,nz./size(z,2),'black');
    xlabel('Z values from the Lorenz attractor')
    ylabel('Probability values')
    
    cd ('../images/')
    filename = strcat('pdf_',num2str(iteration+1));
    print('-djpeg',filename);
    close all
    cd ('../monte_carlo/')
end
matlabpool(close)