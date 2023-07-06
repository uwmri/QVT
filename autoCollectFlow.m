function [] = autoCollectFlow(SVPATH) 
Artery = ["L_ICA";"R_ICA";"L_MCA";"R_MCA";"L_ACA";"R_ACA";"L_PCA";"R_PCA";"BA"];
BranchNumber = [0;0;0;0;0;0;0;0;0];
P = table(Artery,BranchNumber);
fig = uifigure("Position",[100 100 350 375]);
fig.Tag = 'testGUI_tag';
uit = uitable(fig,"Data",P);
uit.ColumnEditable=true;
%ax = uiaxes(fig);
But = uibutton(fig,"Text","Done","Position",[150 340 50 25],...
    "ButtonPushedFcn", @(src,event) plotButtonPushed(fig,uit));
function plotButtonPushed(fig,uit)
writetable(uit.Data,fullfile(SVPATH,'LabelledBranchesQVT.csv'));
close(fig)
end
end