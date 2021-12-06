function refraction
close all; 
% This reads in shallow refraction data in the form of a two column
% file. It defaults to ask for a csv file but space or tab delimited txt or dat files will also work.
%Also xls will work but not the new office 2007 xlsx!!!!!!!!!!
%Column 1 travel time t (ms), Column 2 distance x (m)
%Outputs various useful results!!!
%Version 1.0
%Julian Scott - British Antarctic Survey
%12 February 2008
global error

% Load the default pathname from the saved file (i.e. the last path used)
% Define default pathname as where the mfiles are stored
%thismfilename = which('fbp.txt');     %('fullpath');
%return
% slashind = max(findstr('\',thismfilename));
% this_dir = thismfilename(1:slashind);
% if exist('pathname.mat','file')
%     load pathname
%     if isempty(pathname)
%         pathname = this_dir;
%     end
% else
%     pathname = this_dir;
% end
% 
%     ext = 'csv';
%     if strcmp(pathname(end),'\') || strcmp(pathname(end),'/')
%         pathname = [pathname '*.' ext];
%     else
%         pathname = [pathname '\*.' ext];
%     end
%     
    % Choose file
    pathname = '\Users\epearce\Desktop\fbp.csv';
    this_dir = '\Users\epearce\Desktop\'; 
    [fname,pathname] = uigetfile(pathname,'Choose shallow refraction t,x file');
    % Check if the user cancelled
    if fname == 0 & pathname == 0
        set_message('Loading data cancelled;	Select a data file.');
        return
    else
        % Save the last used path for next time
        full_pathname = [this_dir 'pathname.mat'];
        eval(['save ' full_pathname ' pathname']);
    end
    filename = [pathname fname];


%put data from csv fileahhh into a matrix
%refdata = csvread(filename);

%more generally put data from a 2 column t,x file with no headers into a
%matrix
data = importdata(filename);
refdata(:,1) = data(:,2); %convert to tt,x because that is the way I wrote the program!
refdata(:,2) = data(:,1);

%Sort the matrix by distance (Column 1) in ascending order
[refsort,IX] = sort(refdata(:,1));
refsort(:,2) = refdata(IX,2);


%Now loop through and take average of any offset that has two travel time
%values. Note if the user doesn't want this then don't put the nonrequired travel times
%into the csv file in the first place. A loop may not be the best way of
%doing this but the matrix is fairly small so it is quite quick anyway

m = size(refsort,1);

row = 1;
value = refsort;

for k = 1:m-1
    if refsort(k,1) == refsort(k+1,1)
        value(row,1) = mean(refsort(k:k+1,1));
        value(row,2) = mean(refsort(k:k+1,2));
        row = row+1;
    elseif k>1 & refsort(k,1) == refsort(k-1,1)
        value(row,1) = refsort(k-1,1);
        value(row,2) = refsort(k-1,2);
        row = row+1;
    else
        value(row,1) = refsort(k,1);
        value(row,2) = refsort(k,2);
        row = row + 1;
    end
end

clear refdata

[dummy,rows,columns] = unique(value(:,1));
refdata(rows,2) = value(rows,2);
refdata(rows,1) = value(rows,1);
ind = find(refdata(:,1));
refdata = refdata(ind,:);
x = refdata(:,1);
t = refdata(:,2);

clear refsort refdata ind rows columns dummy value row m

%%% All cleared up and the distance and travel times are now in vectors x and t


%%Define starting values for paramters in a function relating ditance to
%%travel time

% Use the last 20 offsets to give a starting value for Vice = 1/a5
m = size(x,1);
dx = gradient(x);
dt = gradient(t);
V = dx./dt;
Vice = mean(V((m-21):(m-1)));
%ViceError = std(V((m-21):(m-1)));
str = num2str(Vice*1000,'%4.0f');
question = ['From last 20 offsets Vice = ', str,' m/s. This number will go into the minimisation as a starting variable.'...
    ' Would you like to fix Vice instead and not allow the minimization to alter it?'];
button=questdlg(question,...
                'Yes','Yes','No','No');
        if strcmp(button,'Yes')
            prompt = {'Enter velocity in metres per second:'};
            dlg_title = 'Velocity of ice at depth';
            num_lines = 1;
            def = {str};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            Vice = str2double(answer);
            Vice = Vice./1000;
        end


if strcmp(button,'Yes')
    p(1) = 10; p(2) = 0.09; p(3) = 30; p(4) = 0.007; %Use temp parameter p
    tdata = t - x ./Vice; % subtracting parameter a5 because we are not changing it
    model = @analfun1; %analytical function!
    a = fminsearch(model, p);
    a(5) = 1./Vice;
    a = a';
else
    p(1) = 10; p(2) = 0.09; p(3) = 30; p(4) = 0.007; p(5) = 1./Vice; %Use temp parameter p
    tdata = t;
    model = @analfun2; %analytical function!
    a = fminsearch(model, p);
    a = a'
end


clear p

%Give model fitted values at 1 m intervals
xi = 0:x(end);
ti = a(1).*(1 - exp(-a(2).*xi)) + a(3).*(1-exp(-a(4).*xi)) + a(5).*xi

%Now run function zvel to give a velocity depth profile
[v,z] = zvel(a);

%Now run a function to give the density profile
[ro] = zdensity(v,z,Vice);

%Now run a function to calculate two-way travel time
[ztwtt,twtt] = twoway(v,z,Vice);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure plotting section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bdwidth = 40;
topbdwidth = 80;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1  = [bdwidth,... 
	1/2*scnsize(4) + bdwidth,...
	scnsize(3)/2 - 2*bdwidth,...
	scnsize(4)/2 - (topbdwidth + bdwidth)];
pos2 = [pos1(1) + scnsize(3)/2,...
	pos1(2),...
	pos1(3),...
	pos1(4)];
pos3  = [bdwidth,... 
	 bdwidth,...
	scnsize(3)/2 - 2*bdwidth,...
	scnsize(4)/2 - (topbdwidth + bdwidth)];
pos4 = [pos3(1) + scnsize(3)/2,...
	pos3(2),...
	pos3(3),...
	pos3(4)];

figure('Position',pos1)
hold on
plot(x,t,'xk');
plot(xi,ti,'-k');
xlabel('Distance (m)');
ylabel('Travel time (ms)');
hold off

figure('Position',pos2)
hold on
plot(1000.*v,z), grid on
xlabel('Velocity (ms{^{-1}})');
ylabel('Depth (m)');
set(gca,'YDir','reverse')
hold off

figure('Position',pos3)
hold on
plot(ro,z), grid on
xlabel('Density (kgm{^{-3}})');
ylabel('Depth (m)');
set(gca,'YDir','reverse')
hold off

figure('Position',pos4)
hold on
plot(twtt,ztwtt), grid on
xlabel('Two-way travel time (ms)');
ylabel('Depth (m)');
set(gca,'YDir','reverse')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of figure plotting section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save output.mat x t xi ti a v ro z error ztwtt twtt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m1] = size(x,1); [m2] = size(xi,2); [m3] = size(z,2); [m4] = size(ztwtt,2);
dims = [m1 m2 m3 m4];
totalrows = max(dims);
xi=xi'; ti=ti'; v=v'; ro=ro'; z=z'; ztwtt=ztwtt'; twtt=twtt';
if m1 < totalrows
    x(m1+1:totalrows,1) = NaN;
    t(m1+1:totalrows,1) = NaN;
end

if m2 < totalrows
   xi(m2+1:totalrows,1) = NaN;
   ti(m2+1:totalrows,1) = NaN;
end

if m3 < totalrows
   v(m3+1:totalrows,1) = NaN;
   ro(m3+1:totalrows,1) = NaN;
   z(m3+1:totalrows,1) = NaN;
end

if m4 < totalrows
   ztwtt(m4+1:totalrows,1) = NaN;
   twtt(m4+1:totalrows,1) = NaN;
end


OUTPUT(:,1) = x; OUTPUT(:,2) = t;
OUTPUT(:,3) = xi; OUTPUT(:,4) = ti;
OUTPUT(:,5) = z; OUTPUT(:,6) = v.*1000; OUTPUT(:,7) = ro;
OUTPUT(:,8) = ztwtt; OUTPUT(:,9) = twtt;
OUTPUT = OUTPUT';

[file,path] = uiputfile('*.csv','Save Data As');
filename = fullfile(path, file);

fidout = fopen(filename,'w');

    if fidout == -1
       msgbox('Warning: Unable to open the file.')
        error(' ')
    end

       fprintf(fidout, 'A1,A2,A3,A4,A5,\n');
       fprintf(fidout, '%8.5f , %8.5f, %8.5f, %8.5f, %8.5f,\n', a);
       fprintf(fidout, 'Vice = 1/A5 (m/s),\n');
       fprintf(fidout, '%8.0f,\n', Vice.*1000);
       fprintf(fidout, 'Fitting error,\n');
       fprintf(fidout, '%8.4f,\n', error);

       fprintf(fidout, 'Offset (m),tt (ms),xmodel(m),ttmodel (ms),z (m),v (m/s),ro (kgm^(-3),ztwtt (m),twtt (ms),\n');
       fprintf(fidout, '%8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, \n', OUTPUT);

    
fclose(fidout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error,
    function [sse] = analfun1(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        FittedCurve = A.*(1 - exp(-B.*x)) + C.*(1-exp(-D.*x));
        ErrorVector = FittedCurve - tdata;
        sse = sum(ErrorVector .^ 2);
        error = sse;
    end

    function [sse] = analfun2(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        E = params(5);
        FittedCurve = A.*(1 - exp(-B.*x)) + C.*(1-exp(-D.*x)) + E.*x;
        ErrorVector = FittedCurve - tdata;
        sse = sum(ErrorVector .^ 2);
        error = sse;
    end
end


