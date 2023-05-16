function dispersion=lteVarianceCalculation(symbolPUCCH)
%Выделение нужной части ресурсной сетки

%Выделение незанятых пользователями ячеек (элементы, в которых только шум)

%Оценка дисперсии по методу максимального правдопадобия
mathExpectation=sum(symbolPUCCH)/length(symbolPUCCH);
dispersion=abs(sum((symbolPUCCH-mathExpectation).^2))/length(symbolPUCCH);
end