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
% Fama-French 3 Factor European
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


% Risk Model

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


% index of the first day of trading
t0=index(1)
% t0=find(Date==myDate(1))


% Alphas

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
% NumDaysYear=t0-2;

for i=1:T
    if i<=NumDaysYear
        alphamom(i,:)=NaN(1,n);
    else
        ii=i-NumDaysYear;   % start
        iii=i-NumDaysMonth; % end
        momRtn=(tri(iii,:)-tri(ii,:))./tri(ii,:);
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



% Country Constraints
for i=1:n
    myindexlist=allstocks(i).indexlist;
    Num=length(myindexlist);
    
    if Num>1
        for j=1:Num
            if ~isempty(myindexlist(j).index)
                mycountryindex=myindexlist(j).index;
            end
        end
    else
        mycountryindex=myindexlist.index;
    end
    Country{i}=mycountryindex;
end

[countries, ia, ic]=unique(Country);

% number of industries
l=length(countries);

% Boolean matrix
% Every row of matrix F has exactly one entry = 1
% all other entries = 0
F=zeros(n,l);

% Country Dummy
for i=1:n
    F(i,ic(i))=1;
end



% Optimizer

% booksize
booksize=zeros(T,1);

% stock positions
back_weight=zeros(T,n);

% trade
trade=zeros(T,n);

% Daily Profit and Loss (P&L)
PL=zeros(T,1);

% Trade size
tradesize=zeros(T,1);


% replace NaN with 0
tcost(isnan(tcost))=0;
returns(isnan(returns))=0;


mu=0.1;
lambda=0.001;

r_star=300000;
f_star=100000;


for t=t0:T
    CurrentDate=Date(t);
    
    % same covariance matrix each month
    CurrentMonthIndex=find(myDate<=CurrentDate);
    tt=CurrentMonthIndex(end);
    
    w0=back_weight(t-1,:);
    CurrentActive = active(tt).index;
    w=w0(CurrentActive)';
    NumOfStocks=length(CurrentActive); 
    
    Betas = beta(tt).beta';
    COV = covariance(tt).cov;       
    Alpha=alphablend(t-1,CurrentActive)';
    
    Tau=tcost(t-1,CurrentActive)';
    
    myVolume=volume(t-1,CurrentActive)';
    myVolume(isnan(myVolume))=0;
    
    R0=R(CurrentActive',:);
    F0=F(CurrentActive',:);
    
    H=2*mu*[COV -COV; -COV COV];
    g=[2*mu*COV*w-Alpha+lambda*Tau;...
        -2*mu*COV*w+Alpha+lambda*Tau];
    A=[R0' -R0'; -R0' R0'; F0' -F0'; -F0' F0'];
    b=[r_star*ones(m,1)-R0'*w; r_star*ones(m,1)+R0'*w;...
        f_star*ones(l,1)-F0'*w; f_star*ones(l,1)+F0'*w];
    C=[Betas' -Betas'];
    d=-Betas'*w;
    
    LB=zeros(2*NumOfStocks,1);
    
    % Max trade size
    theta=min(myVolume*0.01, 150000);
    
    % Max position size
    pie=min(10*theta, 2.5/100*50000000);
    
    UB=[max(0,min(theta,pie-w));...
        max(0,min(theta,pie+w))];
    
    options = optimset('Algorithm','interior-point-convex');
    options = optimset(options,'Display','iter');
    [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    
    Today=zeros(1,n);
    Yesterday=back_weight(t-1,:);
    
    if (isempty(u))
        trade(t,:)=trade(t-1,:);
    else
        % vector of desired portfolio weights
        y=u(1:NumOfStocks);
        z=u(NumOfStocks+1:end);
        x=y-z+w;        
        Today(CurrentActive)=x;
        trade(t,CurrentActive)=y-z;
    end
    
    trade(t,:)=Today-Yesterday;
    
    back_weight(t,:)=back_weight(t-1,:).*(1+returns(t,:))+trade(t,:);
    
    booksize(t)=sum(abs(back_weight(t,:)));
    
    TransactionCost=tcost(t,:)*abs(trade(t,:))';
    
    PL(t)=back_weight(t-1,:)*returns(t,:)'-TransactionCost;
    
    tradesize(t)=sum(abs(trade(t,:)));
end


% Cumulative Profit and Loss (P&L)
pnl=cumsum(PL);

% annualized Sharpe ratio
R=PL(t0:end)./booksize(t0:end);
AnnualReturn=mean(R)*252;
AnnualVolatility=std(R)*sqrt(252);
sharpe=AnnualReturn/AnnualVolatility


j=1;
% drawdown
for i=t0:T
    High=max(pnl(t0:i));
    if pnl(i)==High
        HighWaterMark(j,1)=High;
        HighWaterMarkIndex(j,1)=i;
        j=j+1;
    end
end

% longest drawdown
% longest time without setting a new high-water mark
longest_dd=max(diff(HighWaterMarkIndex))


% deepest drawdown
% biggest loss from the previous high-water mark
deepest_dd=maxdrawdown(pnl,'arithmetic')




figure(1)
plot(Date(t0:end),pnl(t0:end),'b','linewidth',1)
grid on
title('Cumulative P&L', 'Fontsize',12)
ylabel('Value (€)', 'Fontsize',12)


figure(2)
plot(Date(t0:end),booksize(t0:end),'b','linewidth',1)
grid on
title('Booksize', 'Fontsize',12)
ylabel('Value (€)', 'Fontsize',12)


figure(3)
plot(Date(t0:end),tradesize(t0:end),'b','linewidth',1)
grid on
title('Tradesize', 'Fontsize',12)
ylabel('Value (€)', 'Fontsize',12)


% maximum drawdown plot
[MaxDD, MaxDDIndex] = maxdrawdown(pnl,'arithmetic');

figure(4)
plot(Date(t0:end),pnl(t0:end),'b','linewidth',1)
hold on
plot(Date(MaxDDIndex),pnl(MaxDDIndex),'r-o','MarkerSize',10)
title('Maximum Drawdown', 'Fontsize',12)
ylabel('Cumulative P&L (€)', 'Fontsize',12)
