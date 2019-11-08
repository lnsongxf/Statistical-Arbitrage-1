clc
clear all


% Load data
load database

[T,n]=size(price);

Date=datetime(myday,'InputFormat','dd-MMM-yyyy');  % convert to datetime

% arithmatic return
returns=diff(tri)./tri(1:end-1,:);


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
    myreturns=returns(:,activeindex);
    Rtn=myreturns(select,:);
    
    % replace nan with 0
    Rtn(isnan(Rtn))=0;    
    
    % shrinkage
    [sigma,shrinkage]=cov1para(Rtn);
    shrink(i)=shrinkage;
    covariance(i).cov=sigma;
end



figure(1)
plot(myDate,shrink,'r','linewidth',1)
grid on
ylabel('Shrinkage Intensity', 'Fontsize',14)

