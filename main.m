%% Main

for k=1:1 
    figure(WindowState="maximized",Name=("Image "+ k));
    m=1;

    tempIm=imread("Im"+ k +".bmp");
    
    subplot(2,2,1.5)
    imshow(tempIm)
    title("Image "+ k +" (Size = " + numel(tempIm)/1000 +" KB)")
    
    subplot(2,2,3)
    imshow(rgb2gray(tempIm))
    title("Gray Scale")

    subplot(2,2,4)
    imshow(uint8(double(tempIm(:,:,1))/3+double(tempIm(:,:,2))/3+double(tempIm(:,:,3)/3)))
    title("Avrage Of RGB Channels")

    figure(WindowState="maximized",Name=("Image "+ k));

    for l=[80 50 35 15]
        jp=bmp2jpeg(tempIm,l);
        bm=jpeg2bmp(jp);
        subplot(2,2,m)
        imshow(bm)
        title("Image "+ k +" Quality = "+l)
        xlabel("Size = " + sizeCalc(jp) +"KB")        
        m=m+1;
    end

end




%% Functions


% DownSample And UpSample
function result=downSample2dRate2x2(input) % Function #1  
result=downsample(input,2); % First Dim.
result=(downsample(result',2))'; %2nd Dim.
end


function result=upSample2dRate2x2(input)   % Function #2 
result=upsample(input,2);   %Add Zero
r1=zeros(size(result));
r1(2:size(result,1),:)=result(1:size(result,1)-1,:); % r1= result shift
%                     result  =  result +  r1 
result=result+r1; %   1 1 2 2 = 1 0 2 0 + 0 1 0 2

result=(upsample(result',2)); % 2nd Dim
r1=zeros(size(result));
r1(2:size(result,1),:)=result(1:size(result,1)-1,:);
result=(result+r1)';

end



% Divide Photo Into 8x8 Matrix 

function result = nxmCrop(input,partIndex,n,m)  % Function #3  
d1n=ceil(size(input,1)/n);
a=rem(partIndex-1,d1n);
b=floor((partIndex-1)/d1n);
result= input((n*a)+1:min((n*(a+1)),size(input,1)) ,(m*b)+1:min((m*(b+1)),size(input,2)) );
end

function result = nxmPartNum(inputSize,n,m) % Function #4 
result=ceil(inputSize(1)/n)*ceil(inputSize(2)/m);
end

function result = nxmCropIndex(inputSize,partIndex,n,m) % Function #5   
d1n=ceil((inputSize(1)/n));
a=rem(partIndex-1,d1n);
b=floor((partIndex-1)/d1n);
result=[(n*a)+1,min((n*(a+1)),inputSize(1));(m*b)+1 ,min((m*(b+1)),inputSize(2))];
end

function result= to8x8(input)   % Function #6   
x=16;
m=size(input,1);    n=size(input,2);
a=rem((x-rem(m,x)),x);    b=rem((x-rem(n,x)),x);
result=input;
result(m+1:m+a,:,:)=0;
result(:,n+1:n+b,:)=0;
end


% Quantize and Inverse

function result=quantizeC(inp,qua)  % Function #7   
q=[ 17 18 24 47 99 99 99 99;
    18 21 26 66 99 99 99 99;
    24 26 56 99 99 99 99 99;
    47 66 99 99 99 99 99 99
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99];
r=50.0/qua;
if(qua>50)
    r=2-(qua/50);
end
q=r*q(1:size(inp,1),1:size(inp,2));
result=round(inp./q);
end

function result=iQuantizeC(inp,qua) % Function #8   
q=[ 17 18 24 47 99 99 99 99;
    18 21 26 66 99 99 99 99;
    24 26 56 99 99 99 99 99;
    47 66 99 99 99 99 99 99
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99;
    99 99 99 99 99 99 99 99];

r=50.0/qua;
if(qua>50)
    r=2-(qua/50);
end

q=r*q(1:size(inp,1),1:size(inp,2));
result=(inp.*q);
end

function result=quantizeY(inp,qua)  % Function #9   
q=[ 16 11 10 16 24 40 51 61;
    12 12 14 19 26 58 60 55;
    14 13 16 24 40 57 69 56;
    14 17 22 29 51 87 80 62;
    18 22 37 56 68 109 103 77;
    24 35 55 64 81 104 113 92;
    49 64 78 87 103 121 120 101;
    72 92 95 98 112 100 103 99];
r=50.0/qua;
if(qua>50)
    r=2-(qua/50);
end
q=r*q(1:size(inp,1),1:size(inp,2));
result=round(inp./q);
end

function result=iQuantizeY(inp,qua) % Function #10  
q=[ 16 11 10 16 24 40 51 61;
    12 12 14 19 26 58 60 55;
    14 13 16 24 40 57 69 56;
    14 17 22 29 51 87 80 62;
    18 22 37 56 68 109 103 77;
    24 35 55 64 81 104 113 92;
    49 64 78 87 103 121 120 101;
    72 92 95 98 112 100 103 99];
r=50.0/qua;
if(qua>50)
    r=2-(qua/50);
end
q=r*q(1:size(inp,1),1:size(inp,2));
result=(inp.*q);
end



% ZigZag Scanning And Inverse

function result = zigZagScan(input) % Function #11  
% Init...
a = 1;  b = 1;  n1 = 1;  m1 = 1;    k = 1;
n2 = size(input, 1);  m2 = size(input, 2);
result = zeros(1, n2 * m2);

while ((b <= n2) && (a <= m2))

    if (mod(a + b, 2) == 0)                 % up
        if (b == n1)
            result(k) = input(b, a);        % first line
            if (a == m2)
      	      b = b + 1;
    	    else
                a = a + 1;
            end
            k = k + 1;
        elseif ((a == m2) && (b < n2))   % last column
            result(k) = input(b, a);
            b = b + 1;
            k = k + 1;
        elseif ((b > n1) && (a < m2))    % all other cases
            result(k) = input(b, a);
            b = b - 1;
            a = a + 1;
            k = k + 1;
        end

    else                                    % down
        if ((b == n2) && (a <= m2))       %  last line
            result(k) = input(b, a);
            a = a + 1;
            k = k + 1;

        elseif (a == m1)                   %     first column
            result(k) = input(b, a);
            if (b == n2)
      	      a = a + 1;
    	    else
                b = b + 1;
            end
            k = k + 1;
        elseif ((b < n2) && (a > m1))     % all other cases
            result(k) = input(b, a);
            b = b + 1;
            a = a - 1;
            k = k + 1;
        end
    end
    if ((b == n2) && (a == m2))          % bottom right element
        result(k) = input(b, a);
        break
    end
end
end

function result = iZigZagScan(input, n2, m2)    % Function #12  
% init...
a = 1;  b = 1;  n1 = 1;   m1 = 1;   k = 1;

result = zeros(n2, m2);

while ((b <= n2) && (a <= m2))
    if (mod(a + b, 2) == 0)                % going up
        if (b == n1)
            result(b, a) = input(k);
            if (a == m2)
      	      b = b + 1;
    	    else
                a = a + 1;
            end
            k = k + 1;
        elseif ((a == m2) && (b < n2))
            result(b, a) = input(k);
            b = b + 1;
            k = k + 1;
        elseif ((b > n1) && (a < m2))
            result(b, a) = input(k);
            b = b - 1;
            a = a + 1;
            k = k + 1;
        end

    else                                   % going down
        if ((b == n2) && (a <= m2))
            result(b, a) = input(k);
            a = a + 1;
            k = k + 1;

        elseif (a == m1)
            result(b, a) = input(k);
            if (b == n2)
      	      a = a + 1;
    	    else
                b = b + 1;
            end
            k = k + 1;
        elseif ((b < n2) && (a > m1))
            result(b, a) = input(k);
            b = b + 1;
            a = a - 1;
            k = k + 1;
        end
    end
    if ((b == n2) && (a == m2))
        result(b, a) = input(k);
        break
    end
end
end


% Run-length encoding And Decoding
function result=rle(input)  % Function #13  
num=sum(downsample(input,2));
result=zeros(1,num);
l=1;
c=1;
m=max(size(input,1),size(input,2));
for k=1:m
    if(k==m)
        result(l)=c;
        result(l+1)=input(k);
        l=l+2;
        c=1;
    else
        if(input(k)==input(k+1) && c<64)
            c=c+1;
        else
            result(l)=c;
            result(l+1)=input(k);
            l=l+2;
            c=1;
        end
    end
end
a=ones(1,length(result)/2);
a=upsample(a,2);
result=result-a;
end

function result =rld(input) % Function #14  
k=1; m=1;
n=sum(input(1:2:size(input,2)-1));
result=zeros(1,n);
while k<=size(input,2)
    result(m:m+input(k))=input(k+1);
    m=m+input(k)+1;
    k=k+2;
end
end



% Main Functions


function result= bmp2jpeg(input,quality)    % Function #15  
im1=to8x8(input);
im1ycc=rgb2ycbcr(im1);
im1yccd=double(im1ycc)-128;

im1cbDS=downSample2dRate2x2(im1yccd(:,:,2));
im1crDS=downSample2dRate2x2(im1yccd(:,:,3));

partNum=nxmPartNum(size(im1cbDS),8,8);
cr=zeros(1,partNum*64);
cb=zeros(1,partNum*64);

for k=1:partNum
    tempcb=zigZagScan(quantizeC(dct(nxmCrop(im1cbDS,k,8,8)) , quality));
    tempcr=zigZagScan(quantizeC(dct(nxmCrop(im1crDS,k,8,8)) , quality));

    cb(1+64*(k-1):64*k)= tempcb;
    cr(1+64*(k-1):64*k) = tempcr;
end


im1y=im1yccd(:,:,1);
partNum=nxmPartNum(size(im1y),8,8);
y=zeros(1,partNum*64);

for k=1:partNum
    tempy=zigZagScan(quantizeY(dct(nxmCrop(im1y,k,8,8)) , quality));

    y(1+64*(k-1):64*k)=tempy;
end

base=2048;
m=size(input,1);    n=size(input,2);
m1=floor(m/base);    m2=rem(m,base);  n1=floor(n/base);    n2=rem(n,base);

result=[quality m1 m2 n1 n2 y cb cr];
result=int16(rle(result));
end


function result= jpeg2bmp(input)    % Function #16  
input=double(input);
[quality,mainN,mainM ,y, cb, cr] = yccDecSep(rld(input));

n=mainN+rem((16-rem(mainN,16)),16);    m=mainM+rem((16-rem(mainM,16)),16);
cbIz=zeros(n/2,m/2);
crIz=zeros(n/2,m/2);


for k=1:nxmPartNum([n/2 m/2],8,8) % 8x8 cb and cr
    tempcb=128+idct(iQuantizeC(iZigZagScan(cb(1+64*(k-1):64*k),8,8) , quality));
    tempcr=128+idct(iQuantizeC(iZigZagScan(cr(1+64*(k-1):64*k),8,8) , quality));

    tempInd=nxmCropIndex([n/2 m/2],k,8,8);
    cbIz(tempInd(1):tempInd(3),tempInd(2):tempInd(4))=tempcb;
    crIz(tempInd(1):tempInd(3),tempInd(2):tempInd(4))=tempcr;
end


cbUP=upSample2dRate2x2(cbIz);
crUp=upSample2dRate2x2(crIz);

yIz=zeros(n,m);
for k=1:nxmPartNum([n m],8,8) % 8x8 y
    tempy=128+idct(iQuantizeY(iZigZagScan(y(1+64*(k-1):64*k),8,8) , quality));
    tempInd=nxmCropIndex([n m],k,8,8);
    yIz(tempInd(1):tempInd(3),tempInd(2):tempInd(4))=tempy;
end


result=zeros(mainN,mainM,3);
result(:,:,1)=yIz(1:mainN,1:mainM);
result(:,:,2)=cbUP(1:mainN,1:mainM);
result(:,:,3)=crUp(1:mainN,1:mainM);

result=uint8(ycbcr2rgb(uint8(result)));
end



% Decoding Array Of JPEG 

function [quality,mainN,mainM ,y, cb ,cr]=yccDecSep(input)  % Function #17  
n1=input(2);    n2=input(3);  m1=input(4);    m2=input(5);
base=2048;
quality= input(1);    mainN=base*n1+n2;   mainM=base*m1+m2;
n=mainN+rem((16-rem(mainN,16)),16);    m=mainM+rem((16-rem(mainM,16)),16);

inp=input(6:size(input,2)); %y cb cr elements
s1=ceil(n/2); s2=ceil(m/2); 
a=m*n;  b=s1*s2;    c=a+b;

y=inp(1:a);
cb=inp(a+1:c);
cr=inp(c+1:size(inp,2));
end




% JPEG Size

function result=sizeCalc(input) % Function #18  
    result=round((length(input))*(9/8))/1000;
end



