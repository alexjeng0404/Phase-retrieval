clear; close all; clc;

%% Initialization
N = 2;
Aperture = [1,2;3,4];
Diffraction = fftshift(fft2(Aperture) / N);
Diff_pattern = abs(Diffraction);

% analysis
energy = @(input) trace(input * (input')); 
error = @(input) energy(abs(fftshift(fft2(input)/size(input,1))) - Diff_pattern);

%% Energy Landscape
K = 200;
grid = linspace(0,6,K);
[val1,val2] = meshgrid(grid);
Landscape = zeros(K);

pos1 = [1,1];
pos2 = [2,2];

initial_input = Aperture;
for i = 1:K
    for j = 1:K
        disp(['(i,j) = (',num2str(i),',',num2str(j),')']);

        U = initial_input;
        U(pos1(1),pos1(2)) = val1(i,j);
        U(pos2(1),pos2(2)) = val2(i,j);

        Landscape(i,j) = error(U);
    end
end

%% Plot Energy Landscape
figure;
contour(val1,val2,Landscape,400);
axis square;
set(gcf,'Color',[1,1,1]);
xlabel('$A_{11}$','Interpreter','latex');
ylabel('$A_{22}$','Interpreter','latex');
fontname('Times New Roman');

%% Gradient Descent
figure;
surf(val1,val2,Landscape,'Edgecolor','none');
axis square;
set(gcf,'Color',[1,1,1]);

max_iter = 10;

initial_input = 6 * rand(1,2);
U = [initial_input(1),2;3,initial_input(2)];
hold on
plot3(U(1,1),U(2,2),error(U),'r*','MarkerSize',10)
for i = 1: 20
    Z = fftshift(fft2(U)) / N;
    Y = abs(Z);
    g = 2 * real(ifft2(ifftshift(exp(1i*angle(Z)) .* (Y - Diff_pattern))) * N);
    g(1,2) = 0; g(2,1) = 0;
    
    U = U - 0.1 * g; 
    U = abs(U);
    disp(['error(U) = ', num2str(error(U))]);
    plot3(U(1,1),U(2,2),error(U),'r*','MarkerSize',10)
end
hold off
view(20,40)

xlabel('$A_{11}$','Interpreter','latex');
ylabel('$A_{22}$','Interpreter','latex');
zlabel('Loss');
fontname('Times New Roman');
