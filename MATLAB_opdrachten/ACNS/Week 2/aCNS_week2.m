clear all
% DATA='pixels';
% DATA='gratings';
DATA='natural';

METHOD='pca';
METHOD='sparse';


% Generate data
switch DATA,
case 'pixels'
    n1=16; % image size is n1 * n1
    N=1500; % number of images
    n=n1^2;
    x=zeros(N,n1,n1);
    y=rand(N,n);
    x=-log(y);
    x=reshape(x,N,n1,n1);
    x=x/max(max(max(x)));
case 'gratings'
    n1=16; % image size is n1 * n1
    N=1000; % number of images
    n=n1^2;
    x=zeros(N,n1,n1);
    n2=4;  % number of grating components
    k=2*pi*[1/2,1/4,1/8,1/16]; % grating wave number
    for i=1:N/2,
        y=rand(1,n2);
        a=-log(y);
        x(i,:,:)=ones(n1,1)*(a*sin(k'*(1:n1)));
    end;
    for i=1:N/2,
        y=rand(1,n2);
        a=-log(y);
        x(N/2+i,:,:)=(a*sin(k'*(1:n1)))'*ones(1,n1);
    end;
    x=x/max(max(max(x)));
case 'natural'
    I=imread('sampleMerry_0011_Lasalle.jpeg');
    g=rgb2gray(I);
    figure(5)
    imshow(g)
    [a1,a2]=size(g);
    k=0;
    n1=16;
    n=n1^2;
    N=a1*a2/n;
    x=zeros(N,n1,n1);
    for i=1:a1/n1,
    for j=1:a2/n1,
        k=k+1;
        x(k,:,:)=g((i-1)*n1+1:i*n1,(j-1)*n1+1:j*n1);
    end;
    end;
    x=x/256;
end; 
figure(1)
for k=1:9,
    colormap(gray)
    subplot(3,3,k);
    imagesc(squeeze(x(k,:,:)));
    colorbar;
end;

X=zeros(N,n);
for k=1:N,
    X(k,:)=reshape(squeeze(x(k,:,:)),1,n);
end;
m=mean(X,1);
X=X-ones(N,1)*m;