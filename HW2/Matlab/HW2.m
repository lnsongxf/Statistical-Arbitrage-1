clc
clear all


% Load data
load allstocks
load TickSummary.mat

filename = 'Datastream.xlsm';

% filename = 'cc.xlsm';

[~,sheet_name]=xlsfinfo(filename);

m=1;
for k=2:numel(sheet_name)
    LoadData=readtable(filename,'Sheet',sheet_name{k});
    
    Data{k}=LoadData;
    Num=length(LoadData{4,:})-1;
    N=floor(Num/7);
    for i=1:N
        x=(i-1)*7+1;
        
        % stock name
        stocks(m).name=LoadData{3,x+1};        
        stocksname(m)=LoadData{3,x+1};
        
        % Share price
        SharePrice=LoadData{4:end,x+1};
        stocks(m).price = cellfun(@str2double,SharePrice);
        
        % Total return index
        TotalReturnIndex = LoadData{4:end,x+2};
        stocks(m).tri = cellfun(@str2double,TotalReturnIndex);
        
        % Daily Volume
        DailyValume=LoadData{4:end,x+3};
        stocks(m).volume=cellfun(@str2double,DailyValume);
        
        % Market-to-Book value of equity
        MarketToBook=LoadData{4:end,x+4};
        stocks(m).mtbv=cellfun(@str2double,MarketToBook);
        
        % Total market cap
        MarketCap = LoadData{4:end,x+5};
        stocks(m).cap = cellfun(@str2double,MarketCap);
        
        % Analyst upgrade - Analyst downgrade
        down=LoadData{4:end,x+6};
        up=LoadData{4:end,x+7};
        up = cellfun(@str2double,up);
        down = cellfun(@str2double,down);
        stocks(m).rec=up-down;
        
        m=m+1;
    end
end

myday=datetime(LoadData{4:end,1},'InputFormat','MM/dd/yyyy');  % convert to datetime

stocksname=stocksname';
allnames=unique(stocksname);

% opts = detectImportOptions(filename,'NumHeaderLines',3,'Sheet',sheet_name{2});
% myData = readtable(filename,opts);
% myday=myData{:,1};
% myday.Format = 'dd-MMM-yyy';


% for k=2:numel(sheet_name)
%     opts = detectImportOptions(filename,'NumHeaderLines',3,'Sheet',sheet_name{2});
%     LoadData=readtable(filename,opts);   
%     Data{k}=LoadData;
%     
% end

% myday=LoadData{:,1};
% myday.Format = 'dd-MMM-yyy';

T=length(myday);
n=length(allstocks);

price=zeros(T,n);
tri=zeros(T,n);
volume=zeros(T,n);
mtbv=zeros(T,n);
cap=zeros(T,n);
rec=zeros(T,n);
tcost=zeros(T,n);

for i=1:n
    myindex=allstocks(i).indexlist.index;
end

% myname=allstocks(1).namelist(1).name

id=zeros(1,0);

for i=1:n
    mynamelist=allstocks(i).namelist;
    Number=length(mynamelist);
    
    if Number>1
        for j=1:Number
            if ~isempty(mynamelist(j).name)
                myname=mynamelist(j).name;
            end
        end
    else
        myname=mynamelist.name;
    end
    
%     index=strncmp(myname,stocksname,5);
    index=strcmp(myname,stocksname);
    select=stocks(index);
    
    if ~isempty(select)
        price(:,i)=select.price;
        tri(:,i)=select.tri;
        volume(:,i)=select.volume;
        mtbv(:,i)=select.mtbv;
        cap(:,i)=select.cap;
        rec(:,i)=select.rec;
        id=[id,i];
    else
        price(:,i)=NaN(T,1);
        tri(:,i)=NaN(T,1);
        volume(:,i)=NaN(T,1);
        mtbv(:,i)=NaN(T,1);
        cap(:,i)=NaN(T,1);
        rec(:,i)=NaN(T,1);
    end
end


% average bid-ask spread across all 7 sessions
myspread=mean(spread,3);

myDate=arrayfun(@datestr,datenum(monthlist),'UniformOutput',false);
Date=datetime(myDate);

StartDate=myday(1);
EndDate=myday(end);

% Select Spread Data
ind=find(Date >= StartDate & Date <= EndDate);
SelectDate=Date(ind);
SelectSpread=myspread(ind,:);


% interpolation
NewSpread = interp1(SelectDate,SelectSpread,myday,'previous');

% Transaction Cost
TransactionCost=NewSpread/2;
TransactionCost = fillmissing(TransactionCost,'movmedian',5,'EndValues','nearest');  


for i=1:n
    mydscode=allstocks(i).dscode;
    Number=length(mynamelist);
    index=strcmp(mydscode,dscode);
    
    if sum(index)
        tcost(:,i)=TransactionCost(:,index);
    else
        tcost(:,i)=NaN(T,1);
    end
end


% Construct a boolean matrix
isactivenow=true(T,n);

% for i=1:n
%     myPrice=price(:,i);
%     Number=sum(isnan(myPrice));
%     
%     if Number > 0
%         isactivenow(:,i)=false(T,1);
%     else
%         isactivenow(:,i)=true(T,1);
%     end
% end



for i=1:n
    myPrice=price(:,i);
    Number=sum(isnan(myPrice));
    
    if Number > 0
        isactivenow(:,i)=false(T,1);
    else
        isactivenow(:,i)=true(T,1);
    end
end










save QUEBEC.mat allstocks myday price tri volume mtbv rec tcost isactivenow
