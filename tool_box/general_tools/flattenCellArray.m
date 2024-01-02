% Función recursiva para extraer los elementos no vacíos de un cell array
% Azul Silva, 2023
function flattenedVector = flattenCellArray(cellArray)
    flattenedVector = [];
    for i = 1:numel(cellArray)
        if iscell(cellArray{i})
            flattenedVector = [flattenedVector, flattenCellArray(cellArray{i})];
        else
            if ~isempty(cellArray{i})
                flattenedVector = [flattenedVector, cellArray{i}];
            end
        end
    end
end
