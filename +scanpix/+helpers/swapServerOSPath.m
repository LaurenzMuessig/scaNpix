function pathSwapped = swapServerOSPath(path2swap,options)
%UNTITLED7 Summary of this function goes here

%%
arguments
    path2swap (1,:) {mustBeText}
    options.win (1,:) {mustBeText} = 'S:\';
    options.linux (1,:) {mustBeText} = '/mnt/laurenz/s1/';
end


%%
% swap server path
if strcmp(path2swap(1),'/') && ispc
    pathSwapped = fullfile(replace(path2swap,options.linux,options.win));
elseif isunix && ~strcmp(path2swap(1),'/')
    pathSwapped = fullfile(replace(path2swap,options.win,options.linux));
else
    pathSwapped = path2swap;
end 