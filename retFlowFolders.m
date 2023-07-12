function [Anatpath,APpath,LRpath,SIpath] = retFlowFolders(path2flow,varargin)
    TestDir=dir(path2flow);
    RegExp={'.*FLOW.*(\d{1}\.?\d{1})?.*Anat.*',...
        '.*FLOW.*(\d{1}\.?\d{1})?.*AP.*',...
        '.*FLOW.*(\d{1}\.?\d{1})?.*SI.*',...
        '.*FLOW.*(\d{1}\.?\d{1})?.*LR.*'};
    if ~isempty(varargin)
        res=varargin(1);
        res=res{1};
        splt=strsplit(res,'.');
        if length(splt) == 2
            strcat(splt{1},splt{2});
        end
    else 
        res=[];
        splt=[];
    end
    Anatpath=[];APpath=[];LRpath=[];SIpath=[];
    for i=1:length(TestDir)
        Exp = TestDir(i).name;
        [match] = regexp(Exp,RegExp{1},'tokens');
        if length(match) == 1
            if strcmp(match{1},res) | strcmp(match{1},splt) 
                Anatpath=fullfile(path2flow,Exp);
            elseif isempty(Anatpath)
                Anatpath=fullfile(path2flow,Exp);
            end
        end
        [match] = regexp(Exp,RegExp{2},'tokens');
        if length(match) == 1
            if strcmp(match{1},res) | strcmp(match{1},splt) 
                APpath=fullfile(path2flow,Exp);
            elseif isempty(APpath)
                APpath=fullfile(path2flow,Exp);
            end
        end
        [match] = regexp(Exp,RegExp{3},'tokens');
        if length(match) == 1
            if strcmp(match{1},res) | strcmp(match{1},splt) 
                SIpath=fullfile(path2flow,Exp);
            elseif isempty(SIpath)
                SIpath=fullfile(path2flow,Exp);
            end
        end
        [match] = regexp(Exp,RegExp{4},'tokens');
        if length(match) == 1
            if strcmp(match{1},res) | strcmp(match{1},splt) 
                LRpath=fullfile(path2flow,Exp);
            elseif isempty(LRpath)
                LRpath=fullfile(path2flow,Exp);
            end
        end
    end
end