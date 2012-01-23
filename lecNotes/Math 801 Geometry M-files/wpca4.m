close all
%path = '..\ROOT\GRAVI\ONE ROOT\06.12.27.17.33.00 SCb 3D\';
%d = dir(path);
%N = length(d);
%im = imread(strcat(path, '108.TIF'));
im = imread('108.TIF');
figure(1), imshow(im, []); 
hold on
im0 = im;

F = [1  4  7  4  1
     4 16 26 16  4
     7 26 41 26  7
     4 16 26 16  4
     1  4  7  4  1];
F = F / sum(sum(F));

im1 = im;
for i=1:10
    im1 = filter2(F, double(im1));
end
im3 = im1;
im1(find(im1<255/1.2)) = 0;
im1(find(im1>0)) = 1;

win_width = 10;
win_height = 100;
for i = size(im1, 2) : -1 : win_width+1
    data = im1(:, i-win_width:i);
    [ti tj] = find(data<1); tj = tj + i-win_width-1;
    tj(find(abs(ti-mean(ti))> 2*std(ti))) = [];
    ti(find(abs(ti-mean(ti))> 2*std(ti))) = [];
    ti(find(abs(tj-mean(tj))> 2*std(ti))) = [];
    tj(find(abs(tj-mean(tj))> 2*std(ti))) = [];
    win = [ti tj];
    if size(win, 1) > 300
        win = win - ones(length(ti), 1) * mean(win, 1);
        V = cov(win);
        SD = sqrt(diag(V));
        if (min(SD)>0)
            V = V ./ (SD*SD');
            [COEFF LATENT EXPLAINED] = pcacov(V);
        end
        x1 = mean(tj); y1 = mean(ti);   
        x2 = x1 - 1;     
        if abs(COEFF(1, 2)) > 0
            k = COEFF(2, 2) / COEFF(1, 2);
            y2 = y1 - k*(x1 - x2);
        else
            y2 = y1;
        end
        figure(1), plot([x1 x2], [y1 y2], '--r', 'LineWidth', 1); drawnow;        
    end    
end
hold off
