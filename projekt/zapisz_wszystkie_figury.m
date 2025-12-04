function zapisz_wszystkie_figury()
    % Pobierz wszystkie otwarte figury
    figHandles = findall(0, 'Type', 'figure');
    
    if isempty(figHandles)
        fprintf('Brak otwartych figur do zapisania.\n');
        return;
    end
    
    fprintf('Znaleziono %d figur do zapisania.\n', length(figHandles));
    
    % Dla każdej figury
    for i = 1:length(figHandles)
        fig = figHandles(i);
        
        % Pobierz numer figury
        figNum = fig.Number;
        
        % Pobierz nazwę figury (jeśli jest)
        figName = get(fig, 'Name');
        
        % Stwórz nazwę pliku
        if isempty(figName)
            filename = sprintf('figura_%d.png', figNum);
        else
            % Usuń niedozwolone znaki z nazwy
            figName = regexprep(figName, '[^\w\s-]', '');
            figName = strtrim(figName);
            filename = sprintf('figura_%d_%s.png', figNum, figName);
        end
        
        % Zapisz figurę
        fprintf('Zapisuję: %s\n', filename);
        saveas(fig, filename);
        
        % Alternatywnie, dla lepszej jakości:
        % print(fig, filename, '-dpng', '-r300');  % 300 DPI
    end
    
    fprintf('\nZapisano wszystkie figury!\n');
end