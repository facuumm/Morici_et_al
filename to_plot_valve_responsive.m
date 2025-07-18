%% To plot data as it is in the paper
%% Sorting first shock responssive
t = [-5:0.1:5];
% indices to calculate the mean and sort
[~ , low] = min(abs(t-0));
[~ , up] = min(abs(t-1));

% --- dHPC ---
% Separate the row with IDX == 1 from the rest of the matrix
IDX = curve.id.dHPC(:,2)==1;
row_to_move = curve.dHPC(:,IDX);  % Row where IDX == 1
remaining_rows = curve.dHPC(:,not(IDX));  % The remaining rows where IDX == 0

% Sort both the row with IDX == 1 and the remaining rows by the mean of each row
[~, sort_idx_remaining] = sort(nanmean(remaining_rows(low:up,:)));  % Sort remaining rows by mean
sorted_remaining_rows = remaining_rows(:,sort_idx_remaining);  % Sorted remaining rows

[~, sort_idx_move] = sort(nanmean(row_to_move(low:up,:)));  % Sort row with IDX == 1 by mean (although only one row)
sorted_row_to_move = row_to_move(:,sort_idx_move);  % Sorted row with IDX == 1

% Reconstruct the matrix: first the sorted row with IDX == 1, then the sorted remaining rows
reordered_matrix = [sorted_remaining_rows,sorted_row_to_move];
IDX = logical([zeros(1,length(sort_idx_remaining)),ones(1,length(sort_idx_move))]);

figure
subplot(121)
imagesc(t,[1:size(reordered_matrix,2)],reordered_matrix'),colormap 'gray',clim([0 1])
set(gca, 'YDir','normal');
hold on, xline(0), xline(1)

subplot(122)
hold on;  % Keep the current plot
% Mark the lines (rows in this case) according to the logical vector
for i = 1:length(IDX)
    if IDX(i)  % Check if the logical vector value is true
        plot([1 size(curve.dHPC,2)], [i i], 'r', 'LineWidth', 0.0001);  % Plot a red line on the row
    end
end
ylim([1, length(IDX)])
hold off;  % Release the plot hold

zeros(1,length(sort_idx_remaining))



% --- vHPC ---
% indices to calculate the mean and sort
IDX = curve.id.vHPC(:,2)==1;
row_to_move = curve.vHPC(:,IDX);  % Row where IDX == 1
remaining_rows = curve.vHPC(:,not(IDX));  % The remaining rows where IDX == 0

% Sort both the row with IDX == 1 and the remaining rows by the mean of each row
[~, sort_idx_remaining] = sort(nanmean(remaining_rows(low:up,:)));  % Sort remaining rows by mean
sorted_remaining_rows = remaining_rows(:,sort_idx_remaining);  % Sorted remaining rows

[~, sort_idx_move] = sort(nanmean(row_to_move(low:up,:)));  % Sort row with IDX == 1 by mean (although only one row)
sorted_row_to_move = row_to_move(:,sort_idx_move);  % Sorted row with IDX == 1

% Reconstruct the matrix: first the sorted row with IDX == 1, then the sorted remaining rows
reordered_matrix = [sorted_remaining_rows,sorted_row_to_move];
IDX = logical([zeros(1,length(sort_idx_remaining)),ones(1,length(sort_idx_move))]);

figure
subplot(121)
imagesc(t,[1:size(reordered_matrix,2)],reordered_matrix'),colormap 'gray',clim([0 1])
set(gca, 'YDir','normal');
hold on, xline(0), xline(1)

subplot(122)
hold on;  % Keep the current plot
% Mark the lines (rows in this case) according to the logical vector
for i = 1:length(IDX)
    if IDX(i)  % Check if the logical vector value is true
        plot([1 size(curve.vHPC,2)], [i i], 'r', 'LineWidth', 0.0001);  % Plot a red line on the row
    end
end
ylim([1, length(IDX)])
hold off;  % Release the plot hold

zeros(1,length(sort_idx_remaining))