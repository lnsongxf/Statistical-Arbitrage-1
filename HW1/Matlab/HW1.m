clc
clear all


% Load data
filename = 'Datastream.xlsm';

[~,sheet_name]=xlsfinfo(filename);

for k=2:numel(sheet_name)
    LoadData=readtable(filename,'Sheet',sheet_name{k});
    NameList = LoadData{1,:};
    myData=LoadData(2:end,:);
    myData.Properties.VariableNames=NameList;
    myDate=datetime(myData.DNMC,'InputFormat','MM/dd/yyyy');  % convert to datetime
    myData.DNMC=myDate;
    
    Num=length(myData.Type);
    
    % Add date column    
    if mod(k,5)==2    
        myData.Date=repmat(datetime(1997,12,31),[Num 1]);
    elseif mod(k,5)==3
        myData.Date=repmat(datetime(1998,12,31),[Num 1]);
    elseif mod(k,5)==4
        myData.Date=repmat(datetime(1999,12,31),[Num 1]);
    elseif mod(k,5)==0
        myData.Date=repmat(datetime(2000,12,31),[Num 1]);
    else
        myData.Date=repmat(datetime(2001,12,31),[Num 1]);
    end
    
    % Add index column    
    if k<=6      
        myData.MARKET=repmat({'LBGBEL20'},[Num 1]);
    elseif k<=6+5
        myData.MARKET=repmat({'LCOSE20C'},[Num 1]);
    elseif k<=6+5+5
        myData.MARKET=repmat({'LHEX25IN0801'},[Num 1]);
    elseif k<=6+5+5+5
        myData.MARKET=repmat({'LFSBF120'},[Num 1]);
    elseif k<=6+5+5+5+5
        myData.MARKET=repmat({'LXETRDAX'},[Num 1]);
    elseif k<=6+5+5+5+5+5+5
        myData.MARKET=repmat({'AMSINDX'},[Num 1]);
    elseif k<=6+5+5+5+5+5+5+5+5
        myData.MARKET=repmat({'ITMILAN'},[Num 1]);
    elseif k<=6+5+5+5+5+5+5+5+5+5
        myData.MARKET=repmat({'LOSLOOBX'},[Num 1]);
    elseif k<=6+5+5+5+5+5+5+5+5+5+5
        myData.MARKET=repmat({'LIBEX35I'},[Num 1]);
    elseif k<=6+5+5+5+5+5+5+5+5+5+5+5
        myData.MARKET=repmat({'LSWEDOMX'},[Num 1]);
    else
        myData.MARKET=repmat({'LSMIMIDI'},[Num 1]);
    end
    
    Data{k}=myData;
    
end

M=length(sheet_name);

Table=Data{2};
for i=3:M
    Table=[Table; Data{i}];
end

% Unique dscodes
List=unique(Table.Type);

% length of unique descodes
n=length(List)


for i=1:n
    dscode=List(i);
    index=strcmp(dscode,Table.Type);
    select=Table(index,:);
    
    % dscode
    allstocks(i).dscode=select{1,1};
    
    % namelist
    [mynames, ia]=unique(select.NAME);
    for j=1:length(mynames)
        allstocks(i).namelist(j).name=select.NAME(ia(j));
        allstocks(i).namelist(j).date=select.DNMC(ia(j));
    end
    
    % Industry
    [myindustry, ia]=unique(select.INDXFS);
    for j=1:length(myindustry)
        allstocks(i).industrylist(j).name=select.INDXFS(ia(j));
        allstocks(i).industrylist(j).date=select.Date(ia(j));
    end
    
    % I/B/E/S ticker
    [myIBES, ia]=unique(select.IBTKR);
    for j=1:length(myIBES)
        allstocks(i).ibeslist(j).name=select.IBTKR(ia(j));
        allstocks(i).ibeslist(j).date=select.Date(ia(j));
    end
    
    % Index memberships
    [myIndex, ia]=unique(select.MARKET);
    for j=1:length(myIndex)
        allstocks(i).indexlist(j).name=select.MARKET(ia(j));
        allstocks(i).indexlist(j).date=select.Date(ia(j));
    end    
    
    % ISIN
    [myISIN, ia]=unique(select.ISIN);
    for j=1:length(myISIN)
        allstocks(i).isinlist(j).name=select.ISIN(ia(j));
        allstocks(i).isinlist(j).date=select.Date(ia(j));
    end
    
end


save allstocks.mat allstocks


