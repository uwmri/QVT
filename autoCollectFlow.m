function [] = autoCollectFlow(SVPATH) 
    %SVPATH
    Artery = ["L ICA";"R ICA";"L MCA";"R MCA";"L ACA";"R ACA";"L PCA";"R PCA";"BA"];
    P = table('Size',[9,2],'VariableTypes',{'cellstr','cellstr'},'VariableNames',["Artery","BranchLabel"]);
    P.Artery=Artery;
    if exist(fullfile(SVPATH,'LabelledBranchesQVT.csv')) == 2
        BNums = readLabels(SVPATH);
        P.BranchLabel=BNums(:,3);
    end
    fig = uifigure("Position",[100 100 350 375]);
    fig.Tag = 'testGUI_tag';
    uit = uitable(fig,"Data",P);
    uit.ColumnEditable=true;
    But = uibutton(fig,"Text","Done","Position",[150 340 50 25],...
        "ButtonPushedFcn", @(src,event) plotButtonPushed(fig,uit));
    function plotButtonPushed(fig,uit)
        writetable(uit.Data,fullfile(SVPATH,'LabelledBranchesQVT.csv'));
        close(fig)
    end
end