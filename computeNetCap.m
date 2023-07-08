%% example of how to compute capacity using the SURF non-multiplexed network

% specify data file directory
fileDir = 'SURFNetExp/';
% maximum number of edges that a snapshot can have
kmax = 20;

cap = 0;
totProb = 0;
for k = 0:kmax
    fname = strcat(fileDir,'k',num2str(k),'.txt');
    data = load(fname);
    cap = cap + data(:,1)'*data(:,2);
    totProb = totProb + sum(data(:,2));
end

disp(cap)
