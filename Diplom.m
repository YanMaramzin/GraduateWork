% Данный код предназначен для моделирования PUCCH формата
% и использования алгоритма оценки дисперсии шума в контрольном
% UL канале сети LTE.
clear
close all
NTxAnts = 1;                  % Количество передающих антенн

%Задание параметров для UE
ue = struct;                  % Структура конфигурации UE
ue.NULRB = 6;                 % 6 блоков ресурсов (1,4 МГц)
ue.CyclicPrefixUL = 'Normal'; % Обычный циклический префикс
ue.Hopping = 'Off';           % Отсутствие скачкообразной частоты
ue.NCellID = 150;             % Идентификатор ячейки, указанный в приложении A9 TS36.104
ue.Shortened = 0;             % Отсутствие передачи SRS
ue.NTxAnts = NTxAnts;
ue.NSubframe = 0;
usersPUCCHpower = 3;

SNRdB = [-16.1 10 -8.1 -4.1]; 
SNR = 10^(SNRdB(1)/20);
%% PUCCH 1a настройка
% Бит индикатора гибридного автоматического запроса на повторение (HARQ) установлен равным единице. Только
% для PUCCH 1a требуется один бит
ACK = 1;
pucch = struct;
% Задает размер ресурсов, выделяемых для формата PUCCH 2. Этот параметр
% влияет на процентное расположение передачи PUCCH 1
pucch.ResourceSize = 0;
pucch.DeltaShift = 1;
% Количество циклических сдвигов, используемых для формата PUCCH 1 в блоках ресурсов с
% смесью форматов 1 и 2. Это параметр N1cs, указанный в
% TS36.104, приложение A9
pucch.CyclicShifts = 0;
% Индекс ресурсов для каждого пользователя, которых может быть до 36
usersPUCCHindices =1:36;

%Установка параметров для канала распростаранения
channel = struct;                   
channel.NRxAnts = 1;                % Количество приемных антенн
channel.DelayProfile = 'ETU';       % Профиль задержки канала
channel.DopplerFreq = 70.0;         % Доплеровская частота в Гц
channel.MIMOCorrelation = 'Low';    % Низкая корреляция MIMO
channel.NTerms = 16;                % Генераторы, используемые в модели затухания
channel.ModelType = 'GMEDS';        % Тип модели рэлеевского затухания
channel.InitPhase = 'Random';       % Случайные начальные фазы
channel.NormalizePathGains = 'On';  % Нормализовать мощность профиля задержки
channel.NormalizeTxAnts = 'On';     
channel.InitTime = 0;

% Информация о модуляции SC-FDMA: требуется для получения частоты дискретизации
info = lteSCFDMAInfo(ue);
% Частота дискретизации канала
channel.SamplingRate = info.SamplingRate;   
ueChannelSeed = randi(1000,[1,36]);
% Данный цикл используется для реализации метода Монте-Карло
txACK = 1;

for i=1:36
    % Цикл для моделирования многопользовательского приема
        for user = 1:11
            % Создание пустой сетки ресурсов для UE
            txgrid = lteULResourceGrid(ue);
            % Настройка индекс ресурса для данного UE
            pucch.ResourceIdx = usersPUCCHindices(user);

        % Бит подтверждения для передачи 1-му (целевому) пользователю, пробивной
        % Формат 1 содержит гибридный индикатор ARQ (ARQQ) ACK и для
        % других пользователей он содержит случайный индикатор HARQ. 
        % Поскольку существуетодиночный индикатор, передача будет 
        % осуществляться в формате 1a.DMRS в формате PUCCH 1 не содержит данных.
%         if (mod(user,2)==1)
%             txACK = 1;
%         elseif mod(user,2)==0
%             txACK = 0;
%         end
        % Генерировать PUCCH 1 и его DRS
        % Каждому 
        [pucch1Sym,infoPUCCH] = ltePUCCH1(ue,pucch,txACK);
        pucch1User(:,user)=pucch1Sym;
        pucch1DRSSym =ltePUCCH1DRS(ue,pucch)*10^(usersPUCCHpower/20);
        pucch1Sym=pucch1Sym*10^(usersPUCCHpower/20);


            
        % Создание индексов для PUCCH 1 и DRS, чтобы распределить по
        % сетке ресурсов 
        pucch1Indices = ltePUCCH1Indices(ue,pucch);
        pucch1DRSIndices = ltePUCCH1DRSIndices(ue,pucch);
            
        userOrtSeqIdx(user,:)=infoPUCCH.OrthSeqIdx;
        userAlpha(user,:)=infoPUCCH.Alpha;
        userNCellCyclicShift(user,:)=infoPUCCH.NCellCyclicShift;
        userScrambSeq(user,:)=infoPUCCH.ScrambSeq;
        userSeqGroup(user,:)=infoPUCCH.SeqGroup;
        userNResourceIdx(user,:)=infoPUCCH.NResourceIdx;
        userSeqIdx(user,:)=infoPUCCH.SeqIdx;
        userSeqGroup(user,:)=infoPUCCH.SeqGroup;

        % Отображение PUCCH 1 и PUCCH 1 DRS на сетку ресурсов
        if (~isempty(txACK))
            txgrid(pucch1Indices) = pucch1Sym;
            txgrid(pucch1DRSIndices) = pucch1DRSSym;
        end



        % SC-FDMA модуляция
        txwave = lteSCFDMAModulate(ue,txgrid);

        % Моделирование канала и наложение принятых сигналов.
        % Дополнительные 25 отсчетов, добавленные в конец сигнала.
        % Они предназначены для покрытия диапазона задержек, ожидаемых при
        % моделировании канала (комбинация задержки реализации и
        % разброса задержки канала ). На каждой итерации цикла мы накапливаем
        % сумму каждого переданного сигнала, имитирующую прием
        % всех пользователей на базовой станции.
%         channel.Seed = ueChannelSeed;
%         if (user==1)
%             rxwave = lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);
%         else
%             rxwave = rxwave + lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);
%         end
    end

    %Коэффициент нормировки шума
%     N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0*ue.NTxAnts);
    %Добавление гауссовского шума в приемнике
%     rxwave = awgn(rxwave,100,3);
    rxwave = awgn(txwave,100,3);

    pucch.ResourceIdx = usersPUCCHindices(1);
    offset = lteULFrameOffsetPUCCH1(ue,pucch,rxwave);
        if (offset<25)
            offsetused = offset;
        end

    pucch1Indices = ltePUCCH1Indices(ue,pucch);
    rxgrid2 = lteSCFDMADemodulate(ue,rxwave(:,:));
    

    %Прохождение  CyclicShift
    alpha=userAlpha(1:user,1:4);
    nCS=alpha*12/(2*pi);

[dispersion(i),Pr2(i)]=lteVarianceCalculation(rxgrid2,nCS);
end


user=1:36;
figure
plot(user,mean(dispersion,2)),grid on
xlabel("UE")
ylabel("\sigma^2")
title("Зависимость оценки дисперсии от количества UE")

freeCell=47:-1:12;
Fish=mean(dispersion,2).^2./freeCell;
figure
plot(freeCell,mean(Fish,2)),grid on

% figure
% plot(Pr2),grid on
% figure
% histogram(dispersion(1,:))
% figure
% histfit(dispersion(1,:))
% figure
% histfit(dispersion(2,:))
