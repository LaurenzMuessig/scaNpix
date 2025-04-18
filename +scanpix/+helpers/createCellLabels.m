function labels = createCellLabels(obj)
% create a list of cell labels for a dataset which is useful for some plotting functions.

arguments
    obj (1,:) {mustBeA(obj,'scanpix.ephys')}
end

%%
switch obj.type
    case 'npix'
        labels = strrep(string(strcat('clu_',num2str(obj.cell_ID(:,1)))),' ','');
    case 'dacq'
        labels = string(strcat('c',strip(cellstr(num2str(obj.cell_ID(:,1)))),'t',strip(cellstr(num2str(obj.cell_ID(:,2)))))); 
        endo
    
end