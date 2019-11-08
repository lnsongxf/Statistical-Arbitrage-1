clc
clear all


% Load data
load database

[T,n]=size(price);

Date=datetime(myday,'InputFormat','dd-MMM-yyyy');  % convert to datetime


% arithmatic return
returns=diff(tri)./tri(1:end-1,:);
returns=[zeros(1,n); returns];

% % replace NaN with 0
% returns(isnan(returns))=0;


% (a)
% short-term contrarian
% the mean-reersion alpha with triangular decay

for i=1:n
    myIndustry=allstocks(i).industrylist.industry;
    Industry(i,1)=convertCharsToStrings(myIndustry);
end

[Ind, ia, ic]=unique(Industry);

% number of industries
m=length(Ind);

% Boolean matrix
% Every row of matrix R has exactly one entry = 1
% all other entries = 0
R=zeros(n,m);

% Industry Dummy
for i=1:n
    R(i,ic(i))=1;
end

% Triangular Decay Window
k=0:21;
wts=1/11-1/231*k;
wts=flip(wts(1:end-1));

alpharev=zeros(T,n);

for i=1:T
    if i<=20
        AvgRtn=NaN(1,n);
        alpharev(i,:)=NaN(1,n);
    else
        ii=i-20;
        myRtn=returns(ii:i,:);
        myRtn(isnan(myRtn))=0;
        AvgRtn=-wts*myRtn;
        myAlpha=AvgRtn*(eye(n)-R*inv(R'*R)*R');
        alpharev(i,:)=processing(myAlpha);
    end   
end

w1=0.50;


% (b)
% short-term procyclical
% the analyst recommendations revision alpha
% number of upgrades – number of downgrades

alpharec=zeros(T,n);

% Triangular Decay Window
k=0:45;
wts=1/23-1/1035*k;
wts=flip(wts(1:end-1));

for i=1:T
    if i<=44
        AvgRec=NaN(1,n);
        alpharec(i,:)=NaN(1,n);
    else
        ii=i-44;
        myRec=rec(ii:i,:);
        AvgRec=wts*myRec;
        alpharec(i,:)=processing(AvgRec);
    end   
end

w2=0.25;


% (c)
% long-term contrarian
% value alpha
% book-to-market ratio

btmv=1./mtbv;

alphaval=zeros(T,n);
for i=1:T
    alphaval(i,:)=processing(btmv(i,:));
end

w3=0.15;



% (d)
% Long-term Procylical
% momentum alpha
% straight momentum

alphamom=zeros(T,n);
NumDaysMonth=21;
NumDaysYear=12*NumDaysMonth;

for i=1:T
    if i<=NumDaysYear
        alphamom(i,:)=NaN(1,n);
    else
        ii=i-NumDaysYear;   % start
        iii=i-NumDaysMonth; % end
        momRtn=(tri(iii,:)-tri(ii,:))./tri(ii,:);
%         % replace NaN with 0
%         momRtn(isnan(momRtn))=0;

%         momRtn=prod(returns(ii:iii,:)+1,'omitnan')-1;
        alphamom(i,:)=processing(momRtn);
    end   
end

w4=0.10;


% blend alphas
alphablend=zeros(T,n);

allAlphas=w1*alpharev+w2*alpharec+w3*alphaval+w4*alphamom;

for i=1:T
    alphablend(i,:)=processing(allAlphas(i,:));
end



