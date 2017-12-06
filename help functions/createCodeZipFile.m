function CodeZipFile = createCodeZipFile(zipFilename)
%% create zip file of current code and configuration
fprintf('Creating zip file from current source code in %s\n', zipFilename);

% codeFiles = cell(1,1);
% fid = fopen('init_paths.txt', 'r');
% tline = fgetl(fid);
% i = 1;
% while ischar(tline)
%     codeFiles{i} = tline;
%     tline = fgetl(fid);
%     i = i + 1;
% end
% fclose(fid);

% currentDir = pwd;
% idx = strfind(currentDir, ' Code')-1;
% absPathToCode = [currentDir(1:idx) ' Code'];
absPathToCode = fileparts(pwd);
zip(zipFilename, absPathToCode);

fidZip = fopen(zipFilename);
CodeZipFile = fread(fidZip);
fclose(fidZip);

end