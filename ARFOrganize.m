% Script to sort clumped isotope data from pollycompile.m ARF file.
% Using script by Julia and Andy as starting point.
% Goal is to make this a function for later use in start to finish scripts.

% New additions and changes:
% - Option to filter results by user name
% - Sample averages and SEs saved in structure, not separate file
% - Tweaked some variable names to suit my own brain (sorry...)
% - Now stores d47, d48, D48, d49, D49 values as well.
% - Calculates clumped temps and fluid d18O.

% Want to add/fix:
% - ability to add new ARF file to existing structure, updating older
% entries as appropriate
% - currently a bug when reading in a .txt file
% - Calculate minimum errors using standards

%% Initialize workspace

clear all
close all
fclose all;
%load DataSet.mat

%% Check if new dataset, or updating old structure

% I = input('are you updating an existing dataset (y/n)', '-s');
% I = string(I);
% if strcmp(I,'n') == 1
%     update = 0;
% elseif strcmp(I,'y') == 1
%     update = 1;
% else disp('answer not recognized (y/n)')
% end

%% Read in user selected file. Accepts .xlsx or .txt

[filename,filepath,filterspec]=uigetfile('*.*','Choose an ARF file');
cd(filepath);

if strcmp(filename(end-4:end),'.xlsx') % if you choose an excel file
    [~,~,raw] = xlsread(filename);
    for i=1:length(raw(:,1))
      if strcmp(raw(i,1), 'Samples')
        DataStart = i+1;
        elseif strcmp(raw(i,1),'Untrusted Samples (Flag=0)')
           DataEnd = i-4;
      end
    end
% pull data and format sample names to remove dates and index numbers
    d47 = cell2mat(raw(DataStart:DataEnd,27));
    D47 = cell2mat(raw(DataStart:DataEnd,29));
    d18O = cell2mat(raw(DataStart:DataEnd,25));
    d13C = cell2mat(raw(DataStart:DataEnd,22));
    d48 = cell2mat(raw(DataStart:DataEnd,33));
    D48 = cell2mat(raw(DataStart:DataEnd,35));
    d49 = cell2mat(raw(DataStart:DataEnd,37));
    D49 = cell2mat(raw(DataStart:DataEnd,39));
    USERS = cellstr(raw(DataStart:DataEnd,3));
    RUNID = cellstr(raw(DataStart:DataEnd,2));
    smplchar = char(RUNID);
    SMPL = cellstr(smplchar(:,10:end));

elseif strcmp(filename(end-3:end),'.txt') % if you choose a .txt file
    fid=fopen(filename);
    f = 0;
    col = 200;
    while f==0
        headerrow = fgetl(fid);
        if strcmp(headerrow,'Individual Carbonates')
            headerrow = fgetl(fid); % do this once to get past one row preceeding actual headers
            headers=strread(headerrow,'%s',col,'delimiter','\t');
            f=1;
        end
    end
    raw=textscan(fid,repmat('%s',[1,length(headers)]),'delimiter','\t');
    i = 1;
    while i<length(raw{1})+1
        if strcmp(raw{1}(i),'Samples')
            DataStart=i+1;
        elseif strcmp(raw{1}(i),'Untrusted Samples (Flag=0)')
            DataEnd = i-4;
        end
        i=i+1;
    end
    
    D47 = str2double(raw{29}(DataStart:DataEnd));   %ARF, no acid
    d13C = str2double(raw{22}(DataStart:DataEnd));  %vpdb
    d18O = str2double(raw{25}(DataStart:DataEnd));  %vpdb
    USERS = str2double(raw{3}(DataStart:DataEnd));  %usernames
    RUNID = raw{2}(DataStart:DataEnd); %Sample names
    SMPL = zeros(length(RUNID));
    for j=1:length(RUNID)
        SMPL{j} = RUNID{j}(10:end);
    end
    SMPL = SMPL';
    
else
    disp('file type not recognized (must be .xlsx or .txt')
    return
end

%% clean up sample names for structures below
for i = 1:length(SMPL);
    SMPL{i} = regexprep(SMPL{i},' ',''); % Remove spaces from sample names - KRH 08232017
    SMPL{i} = regexprep(SMPL{i},'\.|-',''); % no dots or dashes in structure names
    if strcmp(SMPL{i}(1), '_') % remove underscore at beginning of names
        SMPL{i} = SMPL{i}(2:end);
    elseif regexp(SMPL{i}(1),'\d') % no numbers in first position, add x
        SMPL{i} = strcat('x',SMPL{i});
    end
end

%% Separate specific user's samples - KRH 08232017
UserName = input('Enter your prepline user name (for everything input ''''):');
USERList = {};
if strcmp(UserName,'') == 1
    USERList = SMPL;
else
    n=1;
    for i = 1:length(SMPL);
        if strcmp(UserName,USERS(i))==1
            USERList{n} = SMPL{i};
            n=n+1;
        end
    end
end

%% create list of unique sample names
SMPLList = unique(char(USERList),'rows');
SMPLList = cellstr(SMPLList);

%% create a structure in which each sample has a sub-structure with isotopic info

% KRH 08232017 - added structure categories for some more values (D48, D49,
% etc)

for i = 1:length(SMPLList);                   
    DataSet.(char(SMPLList(i))).d47reps = NaN;
    DataSet.(char(SMPLList(i))).D47reps = NaN;
    DataSet.(char(SMPLList(i))).d13Creps = NaN;
    DataSet.(char(SMPLList(i))).d18Oreps = NaN;
    DataSet.(char(SMPLList(i))).d48reps = NaN;
    DataSet.(char(SMPLList(i))).D48reps = NaN;
    DataSet.(char(SMPLList(i))).d49reps = NaN;
    DataSet.(char(SMPLList(i))).D49reps = NaN;
    DataSet.(char(SMPLList(i))).d47ave = NaN;
    DataSet.(char(SMPLList(i))).D47ave = NaN;
    DataSet.(char(SMPLList(i))).d13Cave = NaN;
    DataSet.(char(SMPLList(i))).d18Oave = NaN;
    DataSet.(char(SMPLList(i))).d48ave = NaN;
    DataSet.(char(SMPLList(i))).D48ave = NaN;
    DataSet.(char(SMPLList(i))).d49ave = NaN;
    DataSet.(char(SMPLList(i))).D49ave = NaN;
    DataSet.(char(SMPLList(i))).d47se = NaN;
    DataSet.(char(SMPLList(i))).D47se = NaN;
    DataSet.(char(SMPLList(i))).d13Cse = NaN;
    DataSet.(char(SMPLList(i))).d18Ose = NaN;
    DataSet.(char(SMPLList(i))).d48se = NaN;
    DataSet.(char(SMPLList(i))).D48se = NaN;
    DataSet.(char(SMPLList(i))).d49se = NaN;
    DataSet.(char(SMPLList(i))).D49se = NaN;
    DataSet.(char(SMPLList(i))).n = NaN;
    DataSet.(char(SMPLList(i))).T = NaN;
    DataSet.(char(SMPLList(i))).Tse = NaN;
    DataSet.(char(SMPLList(i))).d18Of = NaN;
    DataSet.(char(SMPLList(i))).d18Ofse = NaN;
end

%% sift through samples to match and store properly. Calc Averages and SEs

% KRH 08232017 - added code to process D48, D49, etc.

for j = 1:length(SMPLList);
    for i = 1:length((SMPL));
        if (strcmp((SMPLList(j)),SMPL(i)) == 1); %find sample replicates
            
            %store isotopic data in MyStruct.SampleName.Isotope
            DataSet.(char(SMPLList(j))).D47reps(end+1) = D47(i);
            DataSet.(char(SMPLList(j))).D48reps(end+1) = D48(i);
            DataSet.(char(SMPLList(j))).D49reps(end+1) = D49(i);
            DataSet.(char(SMPLList(j))).d47reps(end+1) = d47(i);
            DataSet.(char(SMPLList(j))).d48reps(end+1) = d48(i);
            DataSet.(char(SMPLList(j))).d49reps(end+1) = d49(i);
            DataSet.(char(SMPLList(j))).d13Creps(end+1) = d13C(i);
            DataSet.(char(SMPLList(j))).d18Oreps(end+1) = d18O(i);
            
            %count how many replicates (values that are not NaN)
            DataSet.(char(SMPLList(j))).n = sum(~isnan(DataSet.(char(SMPLList(j))).D47reps));
            
            %take the average with this new replicate and calc new SE

            DataSet.(char(SMPLList(j))).D47ave = nanmean(DataSet.(char(SMPLList(j))).D47reps);
            DataSet.(char(SMPLList(j))).D47se = nanstd(DataSet.(char(SMPLList(j))).D47reps)/sqrt(DataSet.(char(SMPLList(j))).n);
            
            DataSet.(char(SMPLList(j))).D48ave = nanmean(DataSet.(char(SMPLList(j))).D48reps);
            DataSet.(char(SMPLList(j))).D48se = nanstd(DataSet.(char(SMPLList(j))).D48reps)/sqrt(DataSet.(char(SMPLList(j))).n);

            DataSet.(char(SMPLList(j))).D49ave = nanmean(DataSet.(char(SMPLList(j))).D49reps);
            DataSet.(char(SMPLList(j))).D49se = nanstd(DataSet.(char(SMPLList(j))).D49reps)/sqrt(DataSet.(char(SMPLList(j))).n);

            DataSet.(char(SMPLList(j))).d47ave = nanmean(DataSet.(char(SMPLList(j))).d47reps);
            DataSet.(char(SMPLList(j))).d47se = nanstd(DataSet.(char(SMPLList(j))).d47reps)/sqrt(DataSet.(char(SMPLList(j))).n);
            
            DataSet.(char(SMPLList(j))).d48ave = nanmean(DataSet.(char(SMPLList(j))).d48reps);
            DataSet.(char(SMPLList(j))).d48se = nanstd(DataSet.(char(SMPLList(j))).d48reps)/sqrt(DataSet.(char(SMPLList(j))).n);

            DataSet.(char(SMPLList(j))).d49ave = nanmean(DataSet.(char(SMPLList(j))).d49reps);
            DataSet.(char(SMPLList(j))).d49se = nanstd(DataSet.(char(SMPLList(j))).d49reps)/sqrt(DataSet.(char(SMPLList(j))).n);
            
            DataSet.(char(SMPLList(j))).d13Cave = nanmean(DataSet.(char(SMPLList(j))).d13Creps);
            DataSet.(char(SMPLList(j))).d13Cse = nanstd(DataSet.(char(SMPLList(j))).d13Creps)/sqrt(DataSet.(char(SMPLList(j))).n);
            
            DataSet.(char(SMPLList(j))).d18Oave = nanmean(DataSet.(char(SMPLList(j))).d18Oreps);
            DataSet.(char(SMPLList(j))).d18Ose = nanstd(DataSet.(char(SMPLList(j))).d18Oreps)/sqrt(DataSet.(char(SMPLList(j))).n);
            
        end      
    end    
end

%% Remove NaNs from first rows: - KRH 08232017

for j = 1:length(SMPLList)
    DataSet.(char(SMPLList(j))).D47reps = DataSet.(char(SMPLList(j))).D47reps(2:end);
    DataSet.(char(SMPLList(j))).D48reps = DataSet.(char(SMPLList(j))).D48reps(2:end);
    DataSet.(char(SMPLList(j))).D49reps = DataSet.(char(SMPLList(j))).D49reps(2:end);
    DataSet.(char(SMPLList(j))).d47reps = DataSet.(char(SMPLList(j))).d47reps(2:end);
    DataSet.(char(SMPLList(j))).d48reps = DataSet.(char(SMPLList(j))).d48reps(2:end);
    DataSet.(char(SMPLList(j))).d49reps = DataSet.(char(SMPLList(j))).d49reps(2:end);
    DataSet.(char(SMPLList(j))).d13Creps = DataSet.(char(SMPLList(j))).d13Creps(2:end);
    DataSet.(char(SMPLList(j))).d18Oreps = DataSet.(char(SMPLList(j))).d18Oreps(2:end);
    
    
end
%% Calculate temperatures, fluid d18O values and errors. - KRH 08232017

% These steps use external functions d47toT.m and d18Ofluid.m, which must
% be in working directory.

% Temp calibration set to Kelson et al 2017 (90C), but can be changed using 
% appropriate modifier (see header for d47toT for details).

for i = 1:length(SMPLList)
    [DataSet.(char(SMPLList(i))).T,DataSet.(char(SMPLList(i))).Tse] = D47toT(DataSet.(char(SMPLList(i))).D47ave,'K',DataSet.(char(SMPLList(i))).D47se);
    [DataSet.(char(SMPLList(i))).d18Of,DataSet.(char(SMPLList(i))).d18Ofse] = d18Ofluid(DataSet.(char(SMPLList(i))).d18Oave,DataSet.(char(SMPLList(i))).T,DataSet.(char(SMPLList(i))).Tse);
end

%% Clean up unneeded variables and Save dataset
clear d13C d18O D47 d47 d48 D48 d49 D49 DataEnd DataStart RUNID SMPL smplchar USERS i j n filename filepath filterspec raw
save DataSet.mat

