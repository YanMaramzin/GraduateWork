clear
close all

NTxAnts = 1;                  % Количество передающих антенн

ue = struct;                  % Структура конфигурации UE
ue.NULRB = 6;                 % 6 блоков ресурсов (1,4 МГц)
ue.CyclicPrefixUL = 'Normal'; % Обычный циклический префикс
ue.Hopping = 'Off';           % Отсутствие скачкообразной частоты
ue.NCellID = 150;             % Идентификатор ячейки, указанный в приложении A9 TS36.104
ue.Shortened = 0;             % Отсутствие передачи SRS
ue.NTxAnts = NTxAnts;
ue.NSubframe = 0;
usersPUCCHpower = [3 -10 -3 3];

SNRdB = [-16.1 -12.1 -8.1 -4.1 10]; 
SNR = 10^(SNRdB(5)/20); 
%% PUCCH 1a настройка
% Бит индикатора гибридного автоматического запроса на повторение (HARQ) установлен равным единице. Только
% для PUCCH 1a требуется один бит
ACK = 1;
pucch = struct;  % PUCCH config structure
% Задает размер ресурсов, выделяемых для формата PUCCH 2. Это влияет на
%процентное расположение передачи PUCCH 1
pucch.ResourceSize = 0;
%Параметр PUCCH с дельта-сдвигом, указанный в приложении A 9 TS36.104 [ <#8 1> ]
pucch.DeltaShift = 1;
% Количество циклических сдвигов, используемых для формата PUCCH 1 в блоках ресурсов с
% смесью форматов 1 и 2. Это параметр N1cs, указанный в
% TS36.104, приложение A9
pucch.CyclicShifts = 0;
% Вектор индексов ресурсов PUCCH для всех UES, как указано в TS36.104
% Приложение А9
usersPUCCHindices = [10 30 90 25];

channel = struct;                   % Channel config structure
channel.NRxAnts = 1;                % Количество приемных антенн
channel.DelayProfile = 'ETU';       % Профиль задержки канала
channel.DopplerFreq = 70.0;         % Доплеровская частота в Гц
channel.MIMOCorrelation = 'Low';    % Низкая корреляция MIMO
channel.NTerms = 16;                % Генераторы, используемые в модели затухания
channel.ModelType = 'GMEDS';        % Тип модели рэлеевского затухания
channel.InitPhase = 'Random';       % Случайные начальные фазы
channel.NormalizePathGains = 'On';  % Нормализовать мощность профиля задержки
channel.NormalizeTxAnts = 'On';     % 
channel.InitTime = 0;

% Информация о модуляции SC-FDMA: требуется для получения частоты дискретизации
info = lteSCFDMAInfo(ue);
channel.SamplingRate = info.SamplingRate;   % Channel sampling rate
ueChannelSeed = [10 30 89 5];

for user = 1:2
            % Создать сетку ресурсов
            txgrid = lteULResourceGrid(ue);

            % Configure resource index for this user
            pucch.ResourceIdx = usersPUCCHindices(user);

            %  Бит подтверждения для передачи 1-му (целевому) пользователю, пробивной
            %  Формат 1 содержит гибридный индикатор ARQ (ARQQ) ACK и для
            % других пользователей он содержит случайный индикатор HARQ. Поскольку существует
            %  одиночный индикатор, передача будет осуществляться в формате 1a. 
            % DRS в формате PUCCH 1 не содержит данных.
            if (user==1)
                txACK = ACK;
            else
                txACK = randi([0 1],1,1);
            end
            
            pucch.ResourceIdx = usersPUCCHindices(user);
            % Generate PUCCH 1 and its DRS
            % Different users have different relative power
            [pucch1Sym,infoPUCCH] = ltePUCCH1(ue,pucch,txACK);
            a=abs(pucch1Sym);
            pucch1User(:,user)=pucch1Sym;
            pucch1DRSSym =ltePUCCH1DRS(ue,pucch)* ...
                10^(usersPUCCHpower(user)/20);
            pucch1Sym=pucch1Sym*10^(usersPUCCHpower(user)/20);

            % Generate indices for PUCCH 1 and its DRS
            pucch1Indices = ltePUCCH1Indices(ue,pucch);
            pucch1DRSIndices = ltePUCCH1DRSIndices(ue,pucch);
            
            symbols(user,:)=infoPUCCH.Symbols;
            ortSeqIdx(user,:)=infoPUCCH.OrthSeqIdx;
            alpha(user,:)=infoPUCCH.Alpha;
            NCellCyclicShift(user,:)=infoPUCCH.NCellCyclicShift;
            ScrambSeq(user,:)=infoPUCCH.ScrambSeq;
            SeqGroup(user,:)=infoPUCCH.SeqGroup;

            % Map PUCCH 1 and PUCCH 1 DRS to the resource grid
            if (~isempty(txACK))
                txgrid(pucch1Indices) = pucch1Sym;
                txgrid(pucch1DRSIndices) = pucch1DRSSym;
            end

            % SC-FDMA modulation
            txwave = lteSCFDMAModulate(ue,txgrid);

           % Моделирование канала и наложение принятых сигналов.
           % Дополнительные 25 выборок, добавленные в конец формы сигнала
           %,предназначены для покрытия диапазона задержек, ожидаемых при
           % моделировании канала (комбинация задержки реализации и
           %разброса задержки канала ). На каждой итерации цикла мы накапливаем
           % сумма каждого переданного сигнала, имитирующая прием
           % всеми четырьмя пользователями на базовой станции.
            channel.Seed = ueChannelSeed(user);
            if (user==1)
                rxwave = lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);

            else
                rxwave = rxwave + ...
                    lteFadingChannel(channel,[txwave; zeros(25,NTxAnts)]);
            end
end

N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0*ue.NTxAnts);
noise =2*N*complex(randn(size(rxwave)),randn(size(rxwave)));
rxwave = rxwave + noise;

pucch.ResourceIdx = usersPUCCHindices(1);
offset = lteULFrameOffsetPUCCH1(ue,pucch,rxwave);
        if (offset<25)
            offsetused = offset;
        end


pucch1Indices = ltePUCCH1Indices(ue,pucch);
rxgrid1 = lteSCFDMADemodulate(ue,rxwave(:,:));

figure
surface(abs(rxgrid1)),grid on



ind=ltePUCCH1Indices(ue,pucch);
[re1,reind1] = lteExtractResources(ind,rxgrid1);
txgrid1=lteULResourceGrid(ue);
txgrid1(reind1) = re1;
figure
surface(abs(txgrid1));


pucch.ResourceIdx = usersPUCCHindices(2);
offset = lteULFrameOffsetPUCCH1(ue,pucch,rxwave);
        if (offset<25)
            offsetused = offset;
        end

pucch1Indices = ltePUCCH1Indices(ue,pucch);
rxgrid2 = lteSCFDMADemodulate(ue,rxwave(:,:));
figure
surface(abs(rxgrid2)),grid on
% 
% 
% ind=ltePUCCH1Indices(ue,pucch);
% [re2,reind2] = lteExtractResources(ind,rxgrid1);
% txgrid2=lteULResourceGrid(ue);
% txgrid2(reind2) = re2;
% figure
% surface(abs(txgrid2));




phi=[-1 1 3 -3 3 3 1 1 3 1 -3 3;
      1 1 3 3 3 -1 1 -3 -3 1 -3 3;
      1 1 -3 -3 -3 -1 -3 -3 1 -3 1 -1;
      -1 1 1 1 1 -1 -3 -3 1 -3 3 -1;
      -1 3 1 -1 1 -1 -3 -1 1 -1 1 3;
      1 -3 3 -1 -1 1 1 -1 -1 3 -3 1;
      -1 3 -3 -3 -3 3 1 -1 3 3 -3 1;
      -3 -1 -1 -1 1 -3 3 -1 1 -3 3 1;
      1 -3 3 1 -1 -1 -1 1 1 3 -1 1;
      1 -3 -1 3 3 -1 -3 1 1 1 1 1;
      -1 3 -1 1 1 -3 -3 -1 -3 -3 3 -1;
      3 1 -1 -1 3 3 -3 1 3 1 3 3;
      1 -3 1 1 -3 1 1 1 -3 -3 -3 1;
      3 3 -3 3 -3 1 1 3 -1 -3 3 3;
      -3 1 -1 -3 -1 3 1 3 3 3 -1 1;
       3 -1 1 -3 -1 -1 1 1 3 1 -1 -3;
      1 3 1 -1 1 3 3 3 -1 -1 3 -1;
     -3 1 1 3 -3 3 -3 -3 3 1 3 -1;
     -3 3 1 1 -3 1 -3 -3 -1 -1 1 -3;
     -1 3 1 3 1 -1 -1 3 -3 -1 -3 -1;
     -1 -3 1 1 1 1 3 1 -1 1 -3 -1;
     -1 3 -1 1 -3 -3 -3 -3 -3 1 -1 -3;
      1 1 -3 -3 -3 -3 -1 3 -3 1 -3 3;
      1 1 -1 -3 -1 -3 1 -1 1 3 -1 1;
      1 1 3 1 3 3 -1 1 -1 -3 -3 1;
      1 -3 3 3 1 3 3 1 -3 -1 -1 3;
      1 3 -3 -3 3 -3 1 -1 -1 3 -1 -3;
     -3 -1 -3 -1 -3 3 1 -1 1 3 -3 -3;
     -1 3 -3 3 -1 3 3 -3 3 3 -1 -1;
      3 -3 -3 -1 -1 -3 -1 3 -3 3 1 -1;
];

w=[ 1 1 1 1;
    1 -1 1 -1;
    1 -1 -1 +1;
    ];

baseSeq=exp(1i*phi(1,:)*pi/4);
rAlpha=exp(1i*alpha(1,:));
rAlpha=rAlpha.*transpose(baseSeq);
user1=w(1,:).*rAlpha(:,1:4);
userOne=[];
userOne=[user1(:,1); user1(:,2); user1(:,3); user1(:,4)];

% [corr1,lags]=xcorr(re1(1:48),userOne);
% figure
% plot(lags,corr1)
