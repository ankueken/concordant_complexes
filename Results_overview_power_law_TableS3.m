clear
cd Results_SFAnalysis/
files = dir('hyps*');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subDirs = files(dirFlags);

for d=1:length(subDirs)
    eval(['cd' ' ' subDirs(d).name])

    files_a=dir('analysis_*.txt');
    files_h=dir('hyps_*.txt');

    T_a=table();T_h=table();

    for i=1:length(files_a)
        T_a=[T_a; readtable(files_a(i).name)];
        T_h=[T_h; readtable(files_h(i).name)];
    end

    v = genvarname(['T_' subDirs(d).name])
    eval([v '= [T_a T_h.Strongest T_h.Strong T_h.Weak T_h.Weakest T_h.Super_Weak];']);
    cd ..
end   
    