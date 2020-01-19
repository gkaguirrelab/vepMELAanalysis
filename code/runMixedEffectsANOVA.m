function runMixedEffectsANOVA(responseModality, responseMetric)


%%
dataTable = [];
rowCounter = 1;
for group = 1:length(groups)
    for ss = 1:20
        
        
        dataTable(rowCounter,1) = result.(groups{group}).Melanopsin(ss);
        dataTable(rowCounter,2) = result.(groups{group}).LMS(ss);
        dataTable(rowCounter,3) = result.(groups{group}).LightFlux(ss);
        dataTable(rowCounter,4) = group;
        
        rowCounter = rowCounter + 1;
        
        
    end
end

dataTable = array2table(dataTable);
dataTable.Properties.VariableNames = {'Melanopsin', 'LMS', 'LightFlux', 'Group'};
dataTable.Group = categorical(dataTable.Group);
wsVariable = table([0 1 2]', 'VariableNames', {'Stimulus'});

rm = fitrm(dataTable, 'Melanopsin,LMS,LightFlux~Group', 'WithinDesign', wsVariable);
ranovaTable = ranova(rm);

%% post hoc compare
posthocTable = multcompare(rm,'Stimulus', 'By', 'Group');

end