% jhb_setBetaOdorPaths
try
    %allfolders=dir('C:\\Users\\Jadhavlab\\Documents\\gitRepos\\Jadhav-Lab-Codes\\BetaOdorProject\\figures1-6Functions');
    allfolders=dir('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figures1-6Functions');
    allfolders(1:2)=[];
    allfolders=allfolders([allfolders.isdir]==1);
    allfolders=allfolders(~contains({allfolders.name},'OLD'));

    %addpath('C:\\Users\\Jadhavlab\\Documents\\gitRepos\\Jadhav-Lab-Codes\\BetaOdorProject\\figures1-6Functions');
    addpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figures1-6Functions');
    %rmpath(genpath('C:\\Users\\Jadhavlab\\Documents\\gitRepos\general-repo'));
    rmpath(genpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions\general-repo'));
    for i=1:length(allfolders)
        addpath(genpath(fullfile(allfolders(i).folder,allfolders(i).name)));
    end

catch


    %allfolders=dir('E:\GithubCodeRepositories\Jadhav-Lab-Codes\BetaOdorProject\figures1-6Functions');
    allfolders=dir('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figures1-6Functions');
    allfolders(1:2)=[];
    allfolders=allfolders([allfolders.isdir]==1);
    allfolders=allfolders(~contains({allfolders.name},'OLD'));

    %addpath('E:\GithubCodeRepositories\Jadhav-Lab-Codes\BetaOdorProject\figures1-6Functions');
    addpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figures1-6Functions');
    %rmpath(genpath('E:\\Users\\Jadhavlab\\Documents\\gitRepos\general-repo'));
    rmpath(genpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions\general-repo'));

    for i=1:length(allfolders)
        addpath(genpath(fullfile(allfolders(i).folder,allfolders(i).name)));
    end

end


