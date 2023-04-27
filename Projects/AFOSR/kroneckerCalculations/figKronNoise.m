% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: April 24, 2023

clear; clc; close all;

%% 3D
A = rand(4,4,4);

xgv = 1:n; ygv = 1:n; zgv = 1:n;
[x, y, z] = meshgrid(xgv,ygv,zgv);
scatter3(x(:),y(:),z(:),10,B(:),'filled');
colormap(jet);
colorbar;


%%
n = 20;
R = 20;
A = zeros(n,n,n);
B = zeros(n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            A(i,j,k) = norm([i - n/2, j - n/2, k - n/2; 0 0 0],1);
            B(i,j,k) = norm([i - n/2, j - n/2, k - n/2; 0 0 0],2);
        end
    end
end

% xgv = 1:n; ygv = 1:n; zgv = 1:n;
idxs = 1:n;
[x, y, z] = meshgrid(idxs, idxs, idxs);

figure;
subplot(1,3,1); scatter3(x(:),y(:),z(:),10, B(:),'filled');
colormap(jet); colorbar;

subplot(1,3,2); scatter3(x(:),y(:),z(:),10, A(:),'filled');
% colormap(jet); colorbar;

C = superkron(A,B);
idxs = 1:size(C,1);
[x, y, z] = meshgrid(idxs, idxs, idxs);
subplot(1,3,3); scatter3(x(:),y(:),z(:),10, C(:),'filled');


%%

X = []; Y = []; Z = []; C = [];
for i=1:n
    for j=1:n
        for k=1:n
            X = [X i];
            Y = [Y j];
            Z = [Z k];
            C = [C B(i,j,k)];
        end
    end
end
C = C / max(C);
C = C / std(C);
% C = C * 10;
cc = repmat(C,[3 1])';
figure; scatter3(X,Y,Z,10,cc,'filled'); hold on;
% figure; scatter3(X,Y,Z,1,cc); hold on;

%% Scratch
scatter3(X,Y,Z, [10], [0.5 0.5 0.5])

x = [0.5*X(:); 0.75*X(:); X(:)];
y = [0.5*Y(:); 0.75*Y(:); Y(:)];
z = [0.5*Z(:); 0.75*Z(:); Z(:)];
figure; scatter3(x,y,z)
figure; scatter3(X,Y,Z)

%% Image Example

img = imread('imgIndikaNoise.jpg');
A = img(:,:,1);
figure; imagesc(A)

%% Singal Example

% Image 1
A = sin(0:0.2:2*pi);
A = repmat(A, [size(A,2) ,1]);


% Image 2
B = fspecial('gaussian', size(A), 20); 
% figure; imagesc(B)
% B = A';

C = kron(B,A);
% C = kron(A,B);
figure; 
subplot(1,3,1); imagesc(A); title('A');
subplot(1,3,2); imagesc(B); title('B');
subplot(1,3,3); imagesc(C); title('C');

% imagesc(A)
