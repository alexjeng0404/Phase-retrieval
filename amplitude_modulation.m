clear;
close all;

beta = 0.8;
trial = 100;
img_index = 6;
% 2:smiley face, 4:single-slit
% 5:double-slit, 6:circular
error1 = zeros(1, trial);
error2 = zeros(1, trial);
error3 = zeros(1, trial);
error4 = zeros(1, trial);
error5 = zeros(1, trial);
error6 = zeros(1, trial);

% Read data

% if img_index == 1
%     freq_images = im2double(imread('images.jpg'));
%     freq_images = rgb2gray(freq_images);
% end    

if img_index == 2 % smiley face
    Grid = linspace(-10,10,101);
    [X,Y] = meshgrid(Grid);
    A = ((abs(X)-3).^2 + (Y+3).^2) <= 1;
    A = A + (X.^2 + (Y-1).^2 <= 4^2) .* (Y>1);
    freq_images = fftshift(abs(fft2(A)));
    freq_images = normalize(freq_images);
end

% if img_index == 3 % diamond
%     Grid = linspace(-10,10,101);
%     [X,Y] = meshgrid(Grid);
%     freq_images = abs(X) + abs(Y) < 2;
%     freq_images = double(freq_images);
% end

if img_index == 4 % single-slit
    x = (-2:0.05:2);
    y = (-2:0.05:2);
    A = y.'*x;
    i1 = 0;
    for a = -2:0.05:2
        j1 = 0;
        i1 = i1+1;
        for b = -2:0.05:2
            j1 = j1+1;
            if abs(a) <= 0.2 && abs(b) <= 0.1
                A(i1,j1) = 1;
            else
                A(i1,j1) = 0;
            end
        end
    end
    freq_images = fftshift(abs(fft2(A)));
    freq_images = normalize(freq_images);
end

if img_index == 5 % double-slit
    x = (-2:0.05:2);
    y = (-2:0.05:2);
    A = y.'*x;
    i1 = 0;
    for a = -2:0.05:2
        j1 = 0;
        i1 = i1+1;
        for b = -2:0.05:2
            j1 = j1+1;
            if abs(a) <= 0.2 && abs(b+0.5) <= 0.1
                A(i1,j1) = 1;
            elseif abs(a) <= 0.2 && abs(b-0.5) <= 0.1
                A(i1,j1) = 1;
            else
                A(i1,j1) = 0;
            end
        end
    end
    freq_images = fftshift(abs(fft2(A)));
    freq_images = normalize(freq_images);
end

if img_index == 6 % circle
    Grid = linspace(-10,10,101);
    [X,Y] = meshgrid(Grid);
    A = X.^2 + Y.^2 < 1;
    freq_images = fftshift(abs(fft2(A)));
    freq_images = normalize(freq_images);
end



N = size(freq_images);
a = zeros(N);
a(floor(N(1)/2),floor(N(2)/2)) = 1;
a = normalize(a);
b = a;
c = a;
d = a;
e = a;
f = a;
phi1 = zeros(N);


PA = @(y) proja(y,phi1);
PB = @(y) projb(y, freq_images);

% Gradient Descent
stepsize = 500;
center = [floor((N(1)+1)/2),floor((N(1)+1)/2)];
g = @(input) gradient(freq_images,input);
sqrt_energy = @(field) sqrt(trace(field .* (field'))); 
Normalize = @(input) input * sqrt_energy(freq_images) / sqrt_energy(input);
iter = @(input) Normalize(input - stepsize .* g(input));

% analysis
error = @(input) sqrt_energy(abs(fftshift(fft2(input)/size(input,1))) - freq_images);


% GS loop
for ii = 1:trial

    a = PA(PB(a));
    a = normalize(a);
    error1(ii) = norm((abs(a)) - (freq_images), 'fro');

end


% DM loop
for ii = 1:trial

    b = b + beta * (PB((1+1/beta)*PA(b)-(1/beta)*b) - PA((1-1/beta)*PB(b)+(1/beta)*b));
    b = normalize(b);
    error2(ii) = norm((abs(b)) - (freq_images), 'fro');

end


% RRR loop
for ii = 1:trial

    c = c + 1.8 * (PA(2*PB(c) - c) - PB(c));  % RRR
    c = normalize(c);
    error3(ii) = norm((abs(c)) - (freq_images), 'fro');

end


% revRRR loop
for ii = 1:trial

    d = d + beta * (PB(2*PA(d) - d) - PA(d));
    d = normalize(d);
    error4(ii) = norm((abs(d)) - (freq_images), 'fro');

end


% RAAR loop
for ii = 1:trial

    e = beta*(PB(2*PA(e) - e) + e) + (1-2*beta)*PA(e);
    e = normalize(e);
    error5(ii) = norm((abs(e)) - (freq_images), 'fro');

end


% gradient descent
for ii = 1:trial
    f = iter(f);
    if ii == 1
        f(center(1),center(2)) = 0.25 * (f(center(1)+1,center(2)) + ...
                                         f(center(1)-1,center(2)) + ...
                                         f(center(1),center(2)+1) + ...
                                         f(center(1),center(2)-1));
        Normalize(f);
    end
    f = abs(f);
    Diffraction_result = fftshift(fft2(f));
    Diffraction_result = Normalize(Diffraction_result);
    error6(ii) = norm((abs(Diffraction_result)) - (freq_images), 'fro');
end



% Display results
figure;
subplot(3, 3, 1), imshow(im2gray(abs(A))), title('Aperture');
subplot(3, 3, 2), imshow(im2gray(freq_images)), title('Target Amplitude');
subplot(3, 3, 4), imshow(im2gray(abs(a))), title('Reconstructed by GS');
subplot(3, 3, 5), imshow(im2gray(abs(b))), title('Reconstructed by DM, beta=0.8');
subplot(3, 3, 6), imshow(im2gray(abs(c))), title('Reconstructed by RRR, beta=1.8');
subplot(3, 3, 7), imshow(im2gray(abs(d))), title('Reconstructed by revRRR, beta=0.8');
subplot(3, 3, 8), imshow(im2gray(abs(e))), title('Reconstructed by RAAR, beta=0.8');
subplot(3, 3, 9), imshow(im2gray(abs(Diffraction_result)),[]), title('Reconstructed by gradient descent');


% Plot convergence
figure;
plot(error1, 'k', 'LineWidth', 1.5); hold on;
plot(error2, 'b', 'LineWidth', 1.5);
plot(error3, 'r', 'LineWidth', 1.5);
plot(error4, 'g', 'LineWidth', 1.5);
plot(error5, 'm', 'LineWidth', 1.5);
plot(error6, 'c', 'LineWidth', 1.5);
hold off;

xlabel('Iteration');
ylabel('Relative Error');
legend('GS', 'DM', 'RRR', 'revRRR', 'RAAR','GD');
title('Convergence Comparison');
grid on;

% Functions
function y = proja(y,phi1)
    x = ifft2(ifftshift(y));
    x = abs(x) .* exp(1i * phi1);
    y = fftshift(fft2(x));
    y = normalize(y);
end

function y = projb(y, freq_images)
    y = freq_images .* exp(1i * angle(y));
end

function y = normalize(x)
    y = x / max(abs(x(:)));
end

function g = gradient(target,input)
    N = size(target,1);
    
    Z = fftshift(fft2(input)) / N;
    Y = abs(Z);
   
    g = 2 * real(ifft2(ifftshift(exp(1i*angle(Z)) .* (Y - target))) / N);
end