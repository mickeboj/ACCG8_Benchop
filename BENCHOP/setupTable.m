function input = setupTable(timeVec,errVec)

tBSeuCallUI = timeVec(:,1)
tBSamPutUI = timeVec(:,2)
tBSupoutCallI = timeVec(:,3)
tBSeuCallUII = timeVec(:,4)
tBSamPutUII = timeVec(:,5)
tBSupoutCallII = timeVec(:,6)
rBSeuCallUI = errVec(:,1)
rBSamPutUI = errVec(:,2)
rBSupoutCallI = errVec(:,3)
rBSeuCallUII = errVec(:,4)
rBSamPutUII = errVec(:,5)
rBSupoutCallII = errVec(:,6)
display(tBSeuCallUI)
Methods={'MC','MC-S','QMC-S','MLMC','MLMC-A',...
    'FFT','FGL','COS',...
    'FD','FD-NU','FD-AD',...
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT'};

Table2=table(tBSeuCallUI,tBSamPutUI,tBSupoutCallI,tBSeuCallUII,tBSamPutUII,tBSupoutCallII,'RowNames',Methods)
err=[rBSeuCallUI,rBSamPutUI,rBSupoutCallI,rBSeuCallUII,rBSamPutUII,rBSupoutCallII];
err=round(log10(err));

% Now use this table as input in our input struct:
input.data = Table2;
input.error = err;

% Set the row format of the data values (in this example we want to use
% integers only):
input.dataFormat = {'%.1e'};

% Switch transposing/pivoting your table:
input.transposeTable = 1;

% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 0;

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;

%latex = latexTable(input);
