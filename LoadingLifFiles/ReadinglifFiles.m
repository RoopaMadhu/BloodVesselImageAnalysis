
function [] = ReadinglifFiles(CommonWorkspaceName)
    % This function reads lif files and converts them into .mat workspaces
    % Reading all .lif files in a directory: 
    % This function saves OP in the Movies directory(where all the .lif files are located)by default.
    % If you want to save the OP in the current directory instead, uncomment the line after uigetdir() command.
    
    wd = cd;
%     MoviesDir = cd;
    MoviesDir = uigetdir(); % get the path of directory containing all the .lif files
    %% un comment this line if you want to store the OP of this function in the movies directory 
%     cd(MoviesDir); %move to the Movies Directory
    
    tmp_prompt = questdlg('Do you want analyse only selected lif files in this directory?');
    if strcmp(tmp_prompt,'Yes')
        Files = uigetfile('*.lif','Select the lif you want to analyse','Pick a file','Multiselect','on');
        Files = cellstr(Files);
    else
        dirOP = dir(MoviesDir); %select all the files in the directory
        Files = {dirOP.name};
    end  
%     Files(3)=[];
    
        %% Deleting the non-lif files from the output of dir command
        Fields2Delete = [];
        
        for i = 1:length(Files)
            [~,~,ext] = fileparts(Files{i});
            if  strcmp(ext,'.lif') == 0
                
                Fields2Delete = [Fields2Delete,i];
                
            end
            
        end % end of for loop
            
    Files(Fields2Delete)=[];
    % end of deleting the non-lif files
    
    
    %% Importing the .lif files to matlab
    for i = 1:length(Files)
        
        sprintf('Loading Movie : %d',i);
        
        lif2mat = ci_loadLif(Files{i},'','');
        
        if ~isempty(CommonWorkspaceName)
            save ([CommonWorkspaceName,'_',num2str(i)],'lif2mat','-v7.3');
        else
            tmpname = strsplit(Files{i},'.lif');
            save(tmpname{1},'lif2mat','-v7.3');
        end %of if loop
        
        clear tmpvar tmpname; % clear the tmpvar
        
    end %end of for loop
    
end % of the main function




