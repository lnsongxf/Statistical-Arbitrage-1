clc
clear all


% Load data
load database

[T,n]=size(price);

Date=datetime(myday,'InputFormat','dd-MMM-yyyy');  % convert to datetime

% arithmatic return
returns=diff(tri)./tri(1:end-1,:);
returns=[zeros(1,n); returns];


% Load data
filename = 'Europe_3_Factors_Daily.csv';
opts = detectImportOptions(filename,'NumHeaderLines',6);
FF = readtable(filename,opts);
FF.Properties.VariableNames(1) = {'Date'};

FFDate=datetime(num2str(FF.Date),'InputFormat','yyyyMMdd');  % convert to datetime

% check missing values that are reported as -99.99
checkmissing=find(FF.Mkt_RF==-99.99);
isempty(checkmissing)

Rf=FF.RF/100;
Mkt_Rf=FF.Mkt_RF/100;


% find the first day of the month
% 1997 - 2002
% select > 1998
start=datetime(1998,01,01);
mydays=day(Date);
mymonths=month(Date);
myyears=year(Date);
index=find(diff(mymonths)~=0)+1;

myindex=find(Date(index)>=start);
index=index(myindex);

FirstDays=mydays(index);
myDate=Date(index);
myisactivenow=isactivenow(index,:);

N=length(FirstDays);

shrink=zeros(N,1);

% select active stocks
% use past year of daily returns
for i=1:N
    myrow=index(i);
    select=find((Date>=myDate(i)-years(1)) & (Date<myDate(i)));
    activeindex=find(isactivenow(myrow,:)==1);
    selectdates=Date(select,:);
    myreturns=returns(:,activeindex);
    Rtn=myreturns(select,:);
    
    % replace nan with 0
    Rtn(isnan(Rtn))=0;    
    
    % shrinkage
    [sigma,shrinkage]=cov1para(Rtn);
    shrink(i)=shrinkage;
    covariance(i).cov=sigma;
    active(i).index=activeindex;    
    
    % CAPM Beta
    select2=find((FFDate>=myDate(i)-years(1)) & (FFDate<myDate(i)));
    selectdates2=FFDate(select2,:);
    
    RiskFree=Rf(select2,:);
    Rf_Interp = interp1(selectdates2,RiskFree,selectdates);
    RtnMinusRf=Rtn-Rf_Interp;
    
    MarketExcessRtn=Mkt_Rf(select2,:);
    MktMinusRf = interp1(selectdates2,MarketExcessRtn,selectdates);
    
    NumActive=length(activeindex);
    b=zeros(1,NumActive);
    for j=1:NumActive
        b(j)=regress(RtnMinusRf(:,j),MktMinusRf);
    end
    beta(i).beta=b;
end


