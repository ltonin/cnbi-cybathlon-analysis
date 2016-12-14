function cnbiutil_disp_progress(index, total, pretext)
% cnbiutil_disp_progress(index, total, [pretext])
%
% Given index and total number of iteration, the function display the
% progression as text.

    if nargin < 3
        pretext = '';
    end

    bs = '\b';
    sp = ' ';

    msg = [pretext ' Progress: %d/%d\n'];
    msglen = length(msg) -5;
    if index == 1
    fprintf(1,sp(ones(1,msglen+3)));
    end
    fprintf(1,[bs(mod(0:2*(ceil(log10(index)) + ceil(log10(total)) +msglen)-1,2)+1) msg],index, total);
    
end

