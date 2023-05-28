function [D,Pr2]=lteVarianceCalculation(resourceGrid,nCS)
   %Ортогональные последовательности
    w=[ 1 1 1 1;
        1 -1 1 -1;
        1 -1 -1 1;
        1  1 -1 -1;
        ];
    % Выделение нужной части ресурсной сетки c исключением DMRS
    prb=[resourceGrid(1:12,1:2) resourceGrid(1:12,6:7)];
    prbFFT=fft(prb,12,1);

    % Выделение незанятых пользователями ячеек (элементы, в которых только шум)
    b=mod(int32(nCS)+1,13);
    if size(b,1)>12
        b=b(1:12,:);
    end
    newPRB=[prbFFT(b(:,1),1) prbFFT(b(:,2),2) prbFFT(b(:,3),3) prbFFT(b(:,4),4)];
    if size(b,1)<12
        freeCells=prbFFT(:);
        a=newPRB(:);
        freeCells = freeCells(~ismember(freeCells,a));
        freeCells=reshape(freeCells,size(prbFFT,1)-size(newPRB,1),4);
        newPRB=[newPRB;freeCells];
    end
    
    %Снятие ортогональных последовательностей
    r=newPRB*transpose(w);
    a=abs(r);
    r=r(:);
    r=r(size(nCS,1)+1:length(r));
    
    %Оценка дисперсии по методу максимального правдопадобия
    D=sum(r.*conj(r))/length(r);
    Pr2=-(length(r)/D.^2-3/D.^3*sum(abs(r).^2));
end