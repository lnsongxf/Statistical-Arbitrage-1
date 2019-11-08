function a = processing(m)

% cross-sectionally demean, standardize and windsorize
% alpha = processing(m)
% m(1,n)

n=length(m);
a=zeros(1,n);

% Step 1: Demean
x=m-mean(m,'omitnan');

% Step 2: Standardize
y=x/std(x,'omitnan');

% Step 3: Windsorize
for i=1:n
    if isnan(y(i))
        a(i)=0;
    else
        if abs(y(i))<=3
            a(i)=y(i);
        elseif y(i)>3
            a(i)=3;
        else
            a(i)=-3;
        end
    end
end


