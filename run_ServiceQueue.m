% Script that runs a ServiceQueue simulation many times and plots a
% histogram

%% Set up

% Set up to run 100 samples of the queue.
n_samples = 100;

% Each sample is run up to a maximum time of 1000.
max_time = 480;

% Record how many customers are in the system at the end of each sample.
NInSystemSamples = cell([1, n_samples]);

rng('default')

%% Run the queue simulation

% The statistics seem to come out a little weird if the log interval is too
% short, apparently because the log entries are not independent enough.  So
% the log interval should be long enough for several arrival and departure
% events happen.
for sample_num = 1:n_samples
    q = ServiceQueue(LogInterval=10);
    q.schedule_event(Arrival(1, Customer(1)));
    run_until(q, max_time);
    % Pull out samples of the number of customers in the queue system. Each
    % sample run of the queue results in a column of samples of customer
    % counts, because tables like q.Log allow easy extraction of whole
    % columns like this.
    NInSystemSamples{sample_num} = q.Log.NWaiting + q.Log.NInService;
end

% Join all the samples. "vertcat" is short for "vertical concatenate",
% meaning it joins a bunch of arrays vertically, which in this case results
% in one tall column.
NInSystem = vertcat(NInSystemSamples{:});

% MATLAB-ism: When you pull multiple items from a cell array, the result is
% a "comma-separated list" rather than some kind of array.  Thus, the above
% means
%
%    NInSystem = horzcat(NInSystemSamples{1}, NInSystemSamples{2}, ...)
%
% which horizontally concatenates all the lists of numbers in
% NInSystemSamples.
%
% This is roughly equivalent to "splatting" in Python, which looks like
% f(*args).

%% Make a picture


totallength = length(q.Prob);

R = zeros(1, totallength);
for n = 1:totallength
    R(1, n) = q.Prob{1, n};
end

figr = figure();
tr = tiledlayout(figr, 1, 1);
axr = nexttile(tr);
hold(axr, "on");

Rh = histogram(axr, R, Normalization = "probability", BinMethod = "auto");


% Start with a histogram.  The result is an empirical PDF, that is, the
% area of the bar at horizontal index n is proportional to the fraction of
% samples for which there were n customers in the system.

fig1 = figure();
t1 = tiledlayout(fig1,1,1);
ax1 = nexttile(t1);
h = histogram(ax1, NInSystem, Normalization="probability", BinMethod="integers");
hold(ax1,'on')

theoryprob = [0.527169,0.351446,0.10041,0.01825,0.00243,0.00025];
xvaluesshifted = [0,1,2,3,4,5];
plot(ax1, xvaluesshifted,theoryprob,'o')

% MATLAB-ism: Once you've created a picture, you can use "hold on" to cause
% further plotting function to work with the same picture rather than
% create a new one.
%hold on;

% MATLAB-ism: Once you've created a picture, you can use "hold on" to cause
% further plotting function to work with the same picture rather than
% create a new one.



% For comparison, plot the theoretical results for a M/M/1 queue.
% The agreement isn't all that good unless you run for a long time, say
% max_time = 10,000 units, and LogInterval is large, say 10.
rho = q.ArrivalRate / q.DepartureRate;
P0 = 1 - rho;
nMax = 10;
ns = 0:nMax;
P = zeros([1, nMax+1]);
P(1) = P0;
for n = 1:nMax
    P(1+n) = P0 * rho^n;
end
plot(ns, P, 'o', MarkerEdgeColor='k', MarkerFaceColor='r');


%Generate a histogram for the total time customers spend in the system 
%TotalTime= DepartureTime - ArrivalTime;
%Do we also need to consider RenegeTime-ArrivalTime???
%TotalCustomers=length(q.Served)+length(q.Renegeing)
for i = 1:length(q.Served)
        TotalTimeS(i) = q.Served{i}.DepartureTime - q.Served{i}.ArrivalTime;
end

for i = 1:length(q.Renegeing)
        TotalTimeR(i) = q.Renegeing{i}.RenegeTime - q.Renegeing{i}.ArrivalTime;
end

TotalTime = [TotalTimeS,TotalTimeR];
fig2 = figure();
t2 = tiledlayout(fig2,1,1);
ax2 = nexttile(t2);
time = histogram(ax2,TotalTime, Normalization = 'probability',BinMethod ='auto');

%Generate a hisogram for the time the customers spend waiting in the queue
%BeginService-ArrivalTime
for i = 1:length(q.Served)
        WaitTimeS(i) = q.Served{i}.BeginServiceTime - q.Served{i}.ArrivalTime;
end

for i = 1:length(q.Renegeing)
        WaitTimeR(i) = q.Renegeing{i}.RenegeTime - q.Renegeing{i}.ArrivalTime;
end
WaitTime = [WaitTimeS,WaitTimeR];

fig3 = figure();
t3 = tiledlayout(fig3,1,1);
ax3 = nexttile(t3);
w = histogram(ax3,WaitTime, Normalization= 'probability', BinMethod= 'auto');

%Generate a histogram for the time the customers spend being served 
%DepartureTime-BeginServiceTime 
for i = 1:length(q.Served)
        ServeTime(i) = q.Served{i}.DepartureTime - q.Served{i}.BeginServiceTime;
end
fig4 = figure();
t4 = tiledlayout(fig4,1,1);
ax4 = nexttile(t4);
s = histogram(ax4,ServeTime, Normalization= 'probability', BinMethod= 'auto');



% This sets some paper-related properties of the figure so that you can
% save it as a PDF and it doesn't fill a whole page.
% gcf is "get current figure handle"
% See https://stackoverflow.com/a/18868933/2407278
fig = gcf;
fig.Units = 'inches';
screenposition = fig.Position;
fig.PaperPosition = [0 0 screenposition(3:4)];
fig.PaperSize = [screenposition(3:4)];