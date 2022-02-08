% build all the limbs
numSils = 10;
numPatches  = 72;

limbs = cell(numSils,numPatches);

% loop through silhouettes (image sets)
for iSil = 1:numSils
    
    % loop through points within each silhouette
    for iP = 1:numPatches
        
        limbs{iSil,iP} = [limbRayStarts{iSil}(iP,:); ...
                         limbRayStarts{iSil}(iP+1,:); ...
                         limbRayEnds{iSil}(iP,:); ...
                         limbRayEnds{iSil}(iP+1,:)]';
        
    end
    
end