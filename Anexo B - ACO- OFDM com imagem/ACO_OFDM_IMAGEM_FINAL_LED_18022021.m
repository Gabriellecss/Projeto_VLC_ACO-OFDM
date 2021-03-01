% % Universidade do Estado do Rio de Janeiro - UERJ
% % Faculdade de Engenharia- FEN
% % Engenharia Elétrica com Ênfase em Telecomunicações 
%% PROJETO DE CONCLUSÃO DE CURSO 
%% Comunicação por luz visível (VLC):Análise e síntese por simulação utilizando o MATLAB
% % Gabrielle Cristina de Souza Silva
% % Leonardo Alexandre da Silva
%================================LOG================================
%Última atualização - 01/02/2021
%==============================Descrição============================
%Rotina para obter as curvas de BER x Eb/No do ACO-OFDM com Imagem
%===============================INÍCIO==============================

clear all;                                                                 %LIMPA AS VARIAVEIS DO CODIGO
close all;                                                                 %FECHA TODAS AS FIGURAS
format long;                                                               %AUMENTA O NÚMERO DE CASAS DECIMAIS PARA MAIS DE 15 DIGITOS

% ACO-OFDM

%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%%                                                       +++++   PARAMETROS DE ENTRADA    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

EbNoVec = input('Digite o SNR (1 a 15): ');                                % Entrada de dados (ler do teclado)
M= input('Digite o tipo de QAM (16,32,64,128,256,1024): ');                % Entrada de dados (ler do teclado)
x= input('Digite 0 sem led, 1 com led: ');                                 % Entrada de dados (ler do teclado)
k = log2(M);                                                               %K É O NÚMERO DE BITS POR SIMBOLO
numSymPerFrame = 5000;                                                     %NÚMERO DE QAM SIMBOLOS POR FRAME | !!COMO INTERFERE NA LINHA DO RESHAPE,TEM QUE SER DIVISIVEL POR 8
berEst = zeros(size(EbNoVec));                                             %MATRIX VAZIA PARA O BER ESTIMADO COM 15 LINHAS
cp=16;                                                                     %PREFIXO CICLICO


%******CONSIDERANDO UMA IMAGEM COMO A ENTRADA DO SISTEMA******

im=imread('teladefundo.bmp');                                              %REALIZA A LEITURA DA IMAGEM
im_bin=dec2bin(im(:))';                                                    %TRANSFORMA A INFORMAÇÃO DE DECIMAL PARA BINÁRIO
im_bin=im_bin(:);                                                          %TRANSFORMA TODOS OS ELEMENTOS EM UMA UNICA COLUNA E DIVERSAS LINHAS DEPENDENDO DOS DADOS

sym_rem=mod(k-mod(length(im_bin),k),k);                                    %DEFININDO O MODO DE OPERAÇÃO
padding=repmat('0',sym_rem,1);                                             %RETORNA UM VETOR DE N COPIAS DE A EM LINHA E COLUNA DE TAMANHO B=A*N... B=REPMAT(A,N) 
im_bin_padded=[im_bin;padding];                                            %
cons_data=reshape(im_bin_padded,k,length(im_bin_padded)/k)';               %AJUSTA O TAMANHO DA MATRIZ IM_BIN_PADDED, PELO TAMANHO POR K LINHAS E X COLUNAS
cons_sym_id=bin2dec(cons_data);                                            %CONVERTE OS DADOS DE CONS_DATA DE BINARIO PARA DECIMAL
dataSym=cons_sym_id;                                                       %DATASYM É IGUAL A CONS_SYSM_ID


%_________________________________________________________________________________________________________________________
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
%% VALOR DE BER COM O SNR (dB)

    snrdB = EbNoVec + 10*log10(k)                                          %CONVERTENDO Eb/No to SNR em dB
    numErrs = 0;                                                           %DEFININDO A VARIAVEL DO NUMERO DE ERROS IGUAL A "0"
    numBits = 0;                                                           %DEFININDO A VARIAVEL DO NÚMERO DE BITS IGUAL A "0"
    RR=0;
   
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 1.                                                         +++++ GERANDO QAM    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
     qam_modulated_data =qammod(dataSym,M);                                     %GERANDO QAM (SIMBOLOS,ORDEM), SAÍDA É SÃO NÚMEROS COMPLEXOS
     ytx =qam_modulated_data;                                                   %SAÍDA DO QAM ytx=qam_modulated_data
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 2-3.                                                      +++++   SIMETRIA HERMITIANA    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
     
     %SERIAL PARA PARALELO
     ytxx=ytx.';                                                                 
     
     %HERMITIANA
     y = upsample(ytxx,2);                                                 %INSERÇÃO DE ZEROS NAS SUBPORTADORAS PARES DO ACO-OFDM
     hermitian = [0 y conj(fliplr(y(1:length(y)-1)))];                     %APLICANDO HERMITIANO, ESPELHANDO A MATRIZ Y
     
     %IFFT
     ifftsinal=ifft(hermitian);                                            %APLICANDO A IFFT NO SINAL, SAINDO NÚMEROS REAIS EM COLUNAS NA MATRIZ
     
     %PARALELO PARA SERIAL
     p2sifftsinal=ifftsinal.';

     %ADICIONANDO PREFIXO CICLICO
     cpp2sifftsinal=[p2sifftsinal(end-cp+1:end,:);p2sifftsinal];
     
     if x==0
     
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 8.                                           +++++   REALIZANDO O CORTE DA ONDA (CLIPPING)    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

%      %CLIPPING DO SINAL (CORTANDO VALORES ABAIXO DE ZERO)
     avg=0;                                                                %DEFININDO VARIAVEL AVG=0                                                                       
       clipped=cpp2sifftsinal;                                             %CLIPPED RECEBE OS VALORES DE 'OFDM_SINAL'
       for i=1:length(clipped)                                             %CÓDIGO DE REPETIÇÃO QUE VAI DE 1 ATÉ O COMPRIMENTO DE CLIPPED
 	     if clipped(i) > avg                                               %SE CLIPPED É MAIOR QUE AVG ENTÃO
 		 clipped(i) = clipped(i);                                          %CLIPPED RECEBE O VALOR DO PROPRIO CLIPPED
     
         elseif clipped(i) < -avg                                          %SE CLIPPED < 0 ENTÃO
 		    clipped(i) = 0;                                                %CLIPPED RECEBE O VALOR "0"
         end                                                               %ENCERRA O COMANDO END
       end                                                                 %ENCERRA O COMANDO FOR

       
     %SINAL RECEBIDO+RUÍDO AWGN   
     sinalrecebido = awgn(clipped,snrdB,'measured');                       %SINAL RECEBIDO SENDO APLICADO O RUÍDO AWGN
       
     %REMOVENDO PREFIXO CICLICO
     rmcpsinalrecebido=sinalrecebido(cp+1:end,:);                          %SINAL RECEBIDO SEM O PREFIXO CICLICO

     
     %SERIAL PARA PARALELO

     s2psinalrecebido=rmcpsinalrecebido.';                                 %MODIFICANDO O DADO PARA PARALELO
     
     %FFT
     
     fftsinal=fft(s2psinalrecebido);                                       %SINAL COM FFT
     
     %PARALELO PARA SERIAL
     
     p2sfftsinal=fftsinal.';                                               %DADO DE PARALELO PARA SERIAL
     
     %SIMETRIA HERMITIANA DEMOD

     SMHp2sfftsinal = p2sfftsinal(2:(length(p2sfftsinal)-1)/2+01);                             
     SMHp2sfftsinal_final = 2*downsample(SMHp2sfftsinal,2);
     
     
     %DEMODULAÇÃO
     qamdemodulado = qamdemod(SMHp2sfftsinal_final,M);                     %DEMODULA O SINAL RECEBIDO, EM SERIAL, NÚMEROS DECIMAIS
     dataOut = de2bi(qamdemodulado,'left-msb');                            %DADO RECEBIDO EM BINÁRIO
     dataout2=fliplr(dataOut);                                             %ESPELHANDO A MATRIZ RECEBIDA
     numBits = numBits + numSymPerFrame*k;                                 %INCREMENTO DA QUANTIDADE DE BITS ENVIADOS
     
     
% ******* Plot the simulation result *********%
f1 = figure(1);                                                       %DEFININDO COMO FIGURA 12
    
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% XX.                                                        +++++     RECEPÇÃO DA IMAGEM TRANSMITIDA ++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

%PASSANDO PARA BINARIO E REMOVENDO SIMBOLO PADDING
rec_syms_cons=dec2bin(qamdemodulado);                               %CONVERSÃO DE DECIMAL PARA BINARIO
rec_im_bin=reshape(rec_syms_cons',numel(rec_syms_cons),1);                 %AJUSTA O TAMANHO DE REC_IM_BIN DE ACORDO COM AS VARIAVEIS, NUMEL PEGA O TAMANHO DA MATRIZ
rec_im_bin=rec_im_bin(1:end-sym_rem);                                      %
ber=sum(abs(rec_im_bin-im_bin))/length(im_bin);                            %CALCULA O VALOR DO BER QUE É IGUAL AO VALOR RECEBIDO MENOS O ORIGINAL SOBRE O ORIGINAL


%% RECUPERANDO A IMAGEM
rec_im=reshape(rec_im_bin,8,numel(rec_im_bin)/8);                          %AJUSTA DO TAMANHO DE REC_IM_BIN EM RELAÇÃO A LINHA E COLUNA, SENDO QUE NUMEL RETORNA O NÚMERO DE ELEMENTOS
rec_im=uint8(bin2dec(rec_im'));                                            %TRANSFORMA DE BINARIO PARA DECIMAL E ALÉM DISSO CONVERTE OS VALORES DE REC_IM NO TIPO UINT8
rec_im=reshape(rec_im,size(im));                                           %AJUSTA O TAMANHO DA MATRIZ REC_IM PARA O TAMANHO DE "IM"


%% GERANDO PLOTAGENS

%CONSTELAÇÃO TRANSMITIDA



% *******GRÁFICO FIGURA 9 - PLOTAGEM GRÁFICO DA CONSTELAÇÃO*******

n_taps=8;                                                                  %
subplot(2,2,1);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS, 2 COLUNAS E O GRÁFICO NA POSIÇÃO 1
grid on                                                                    %ATIVANDO LINHA DE GRADE
plot(qam_modulated_data,'o','linewidth',2,'markersize',10);                %PLOTAGEM DO GRÁFICO DO SINAL MODULADO
xlabel('Fase');                                                            %RÓTULO DO EIXO X DO GRÁFICO
ylabel('Quadratura');                                                      %RÓTULO DO EIXO Y DO GRÁFICO
title('MODULAÇÃO TRANSMITIDA');                                            %TÍTULO DO GRÁFICO
xlim([-10 10]);                                                            %LIMITE DO TAMANHO DO EIXO X DO GRÁFICO
ylim([-10 10]);                                                            %LIMITE DO TAMANHO DO EIXO Y DO GRÁFICO
title(sprintf('\\bfCONSTELAÇÃO TRANSMITIDA\n\\rm%s Modulação: %d',M));       %TÍTULO DO GRÁFICO PLOTADO
grid on;                                                                   %LINHA DE GRADE DO GRÁFICO

%CONSTELAÇÃO RECUPERADA

% *******GRÁFICO FIGURA 9 - PLOTAGEM GRÁFICO DA CONSTELAÇÃO RECUPERADA*******

subplot(2,2,2);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS, 2 COLUNAS E O GRÁFICO NA POSIÇÃO 2
plot(qamdemodulado,'o','markersize',7);                        %PLOTAGEM DO GRÁFICO DO SINAL RECEBIDO
%plot(qamdemodulado(1:500:end),'0','markersize',3);                        %PLOTAGEM DO GRÁFICO DO SINAL RECEBIDO
xlabel('Fase');                                                            %RÓTULO DO EIXO X
ylabel('Quadratura');                                                      %RÓTULO DO EIXO Y
xlim([-10 1000]);                                                            %LIMITAÇÃO DO EIXO X DO GRÁFICO
ylim([-1 15]);                                                            %LIMITAÇÃO DO EIXO Y DO GRÁFICO
title(sprintf('\\bfCONSTELAÇÃO RECEBIDA\n\\rmVALOR DO SNR: %.2f dB: %s',snrdB,k)); %TÍTULO DO GRÁFICO DO GRÁFICO DA CONSTELAÇÃO RECUPERADA
grid on;                                                                   %LINHA DE GRADE DO GRÁFICO


%IMAGEM ORIGINAL

% *******FIGURA 9 - PLOTAGEM DA IMAGEM ORIGINAL*******

subplot(2,2,3);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS E 2 COLUNAS E O GRÁFICO NA POSIÇÃO 3
imshow(im);                                                                %EXIBE A IMAGEM ORIGINAL
title(sprintf('\\bfIMAGEM TRANSMITIDA\n\\rm: %.2g'));                      %TITULO DO GRÁFICO CRIADO

%IMAGEM RECUPERADA

% *******FIGURA 9 - PLOTAGEM DA IMAGEM RECUPERADA*******

subplot(2,2,4);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS E 2 COLUNAS E O GRÁFICO NA POSIÇÃO 4
imshow(rec_im);                                                            %EXIBE A IMAGEM RECUPERADA NA RECEPÇÃO DESTE PROCESSO                                     
title(sprintf('\\bfIMAGEM RECUPERADA\n\\rmBER: %.2g',ber));                %TÍTULO DO GRÁFICO PLOTADO

%CONFIGURAÇÕES EXTRAS RELACIONADAS A IMAGEM
set(gcf,'Position',[680 287 698 691]);                                     %CONFIGURAÇÃO DO POSICIONAMENTO DA FIGURA
save_file=0;                                                               %OPÇÃO PARA SALVAR O ARQUIVO

sPlotFig = scatterplot(SMHp2sfftsinal_final,1,0,'g.');                     %PLOTAGEM DO SINAL RECEBIDO COM RUÍDO
hold on                                                                    %RETENDO A PLOTAGEM FEITA
scatterplot(ytx,1,0,'k*',sPlotFig)                                         %PONTOS MQAM PLOTADOS NO GRÁFICO EM SOBREPOSIÇÃO AO DADO RECEBIDO


     else
         
     
     s=cpp2sifftsinal;
     p=1;
     imax=1;
     Vt = 0.5;
     Is = 0.5;
     iLED = Is*(exp(s/Vt)-1);                                              %WILLIAM SHOCKLEY`S DIODE EQUATION, corrente do LED
     h=((1+(iLED/imax).^(2*p)).^(1/(2*p)));
     xLED = iLED./h;
       
     
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 8.                                           +++++   REALIZANDO O CORTE DA ONDA (CLIPPING)    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

%      %CLIPPING DO SINAL (CORTANDO VALORES ABAIXO DE ZERO)
     avg=0;                                                                %DEFININDO VARIAVEL AVG=0                                                                       
       clipped=xLED;                                             %CLIPPED RECEBE OS VALORES DE 'OFDM_SINAL'
       for i=1:length(clipped)                                             %CÓDIGO DE REPETIÇÃO QUE VAI DE 1 ATÉ O COMPRIMENTO DE CLIPPED
 	     if clipped(i) > avg                                               %SE CLIPPED É MAIOR QUE AVG ENTÃO
 		 clipped(i) = clipped(i);                                          %CLIPPED RECEBE O VALOR DO PROPRIO CLIPPED
     
         elseif clipped(i) < -avg                                          %SE CLIPPED < 0 ENTÃO
 		    clipped(i) = 0;                                                %CLIPPED RECEBE O VALOR "0"
         end                                                               %ENCERRA O COMANDO END
       end                                                                 %ENCERRA O COMANDO FOR

       
     %SINAL RECEBIDO+RUÍDO AWGN   
     sinalrecebido = awgn(clipped,snrdB,'measured');                       %SINAL RECEBIDO SENDO APLICADO O RUÍDO AWGN
       
     %REMOVENDO PREFIXO CICLICO
     rmcpsinalrecebido=sinalrecebido(cp+1:end,:);                          %SINAL RECEBIDO SEM O PREFIXO CICLICO

     
     %SERIAL PARA PARALELO

     s2psinalrecebido=rmcpsinalrecebido.';                                 %MODIFICANDO O DADO PARA PARALELO
     
     %FFT
     
     fftsinal=fft(s2psinalrecebido);                                       %SINAL COM FFT
     
     %PARALELO PARA SERIAL
     
     p2sfftsinal=fftsinal.';                                               %DADO DE PARALELO PARA SERIAL
     
     %SIMETRIA HERMITIANA DEMOD

     SMHp2sfftsinal = p2sfftsinal(2:(length(p2sfftsinal)-1)/2+01);                             
     SMHp2sfftsinal_final = 2*downsample(SMHp2sfftsinal,2);
     
     
     %DEMODULAÇÃO
     qamdemodulado = qamdemod(SMHp2sfftsinal_final,M);                     %DEMODULA O SINAL RECEBIDO, EM SERIAL, NÚMEROS DECIMAIS
     dataOut = de2bi(qamdemodulado,'left-msb');                            %DADO RECEBIDO EM BINÁRIO
     dataout2=fliplr(dataOut);                                             %ESPELHANDO A MATRIZ RECEBIDA
     numBits = numBits + numSymPerFrame*k;                                 %INCREMENTO DA QUANTIDADE DE BITS ENVIADOS
     
     
% ******* Plot the simulation result *********%
f1 = figure(1);                                                       %DEFININDO COMO FIGURA 12
    
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% XX.                                                        +++++     RECEPÇÃO DA IMAGEM TRANSMITIDA ++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

%PASSANDO PARA BINARIO E REMOVENDO SIMBOLO PADDING
rec_syms_cons=dec2bin(qamdemodulado);                               %CONVERSÃO DE DECIMAL PARA BINARIO
rec_im_bin=reshape(rec_syms_cons',numel(rec_syms_cons),1);                 %AJUSTA O TAMANHO DE REC_IM_BIN DE ACORDO COM AS VARIAVEIS, NUMEL PEGA O TAMANHO DA MATRIZ
rec_im_bin=rec_im_bin(1:end-sym_rem);                                      %
ber=sum(abs(rec_im_bin-im_bin))/length(im_bin);                            %CALCULA O VALOR DO BER QUE É IGUAL AO VALOR RECEBIDO MENOS O ORIGINAL SOBRE O ORIGINAL


%% RECUPERANDO A IMAGEM
rec_im=reshape(rec_im_bin,8,numel(rec_im_bin)/8);                          %AJUSTA DO TAMANHO DE REC_IM_BIN EM RELAÇÃO A LINHA E COLUNA, SENDO QUE NUMEL RETORNA O NÚMERO DE ELEMENTOS
rec_im=uint8(bin2dec(rec_im'));                                            %TRANSFORMA DE BINARIO PARA DECIMAL E ALÉM DISSO CONVERTE OS VALORES DE REC_IM NO TIPO UINT8
rec_im=reshape(rec_im,size(im));                                           %AJUSTA O TAMANHO DA MATRIZ REC_IM PARA O TAMANHO DE "IM"


%% GERANDO PLOTAGENS

%CONSTELAÇÃO TRANSMITIDA



% *******GRÁFICO FIGURA 9 - PLOTAGEM GRÁFICO DA CONSTELAÇÃO*******

n_taps=8;                                                                  %
subplot(2,2,1);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS, 2 COLUNAS E O GRÁFICO NA POSIÇÃO 1
grid on                                                                    %ATIVANDO LINHA DE GRADE
plot(qam_modulated_data,'o','linewidth',2,'markersize',10);                %PLOTAGEM DO GRÁFICO DO SINAL MODULADO
xlabel('Fase');                                                            %RÓTULO DO EIXO X DO GRÁFICO
ylabel('Quadratura');                                                      %RÓTULO DO EIXO Y DO GRÁFICO
title('MODULAÇÃO TRANSMITIDA');                                            %TÍTULO DO GRÁFICO
xlim([-10 10]);                                                            %LIMITE DO TAMANHO DO EIXO X DO GRÁFICO
ylim([-10 10]);                                                            %LIMITE DO TAMANHO DO EIXO Y DO GRÁFICO
title(sprintf('\\bfCONSTELAÇÃO TRANSMITIDA\n\\rm%s Modulação: %d',M));       %TÍTULO DO GRÁFICO PLOTADO
grid on;                                                                   %LINHA DE GRADE DO GRÁFICO

%CONSTELAÇÃO RECUPERADA

% *******GRÁFICO FIGURA 9 - PLOTAGEM GRÁFICO DA CONSTELAÇÃO RECUPERADA*******

subplot(2,2,2);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS, 2 COLUNAS E O GRÁFICO NA POSIÇÃO 2
plot(qamdemodulado,'o','markersize',7);                        %PLOTAGEM DO GRÁFICO DO SINAL RECEBIDO
%plot(qamdemodulado(1:500:end),'0','markersize',3);                        %PLOTAGEM DO GRÁFICO DO SINAL RECEBIDO
xlabel('Fase');                                                            %RÓTULO DO EIXO X
ylabel('Quadratura');                                                      %RÓTULO DO EIXO Y
xlim([-10 1000]);                                                            %LIMITAÇÃO DO EIXO X DO GRÁFICO
ylim([-1 15]);                                                            %LIMITAÇÃO DO EIXO Y DO GRÁFICO
title(sprintf('\\bfCONSTELAÇÃO RECEBIDA\n\\rmVALOR DO SNR: %.2f dB: %s',snrdB,k)); %TÍTULO DO GRÁFICO DO GRÁFICO DA CONSTELAÇÃO RECUPERADA
grid on;                                                                   %LINHA DE GRADE DO GRÁFICO


%IMAGEM ORIGINAL

% *******FIGURA 9 - PLOTAGEM DA IMAGEM ORIGINAL*******

subplot(2,2,3);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS E 2 COLUNAS E O GRÁFICO NA POSIÇÃO 3
imshow(im);                                                                %EXIBE A IMAGEM ORIGINAL
title(sprintf('\\bfIMAGEM TRANSMITIDA\n\\rm: %.2g'));                      %TITULO DO GRÁFICO CRIADO

%IMAGEM RECUPERADA

% *******FIGURA 9 - PLOTAGEM DA IMAGEM RECUPERADA*******

subplot(2,2,4);                                                            %REPRESENTA UM QUADRANTE DE 2 LINHAS E 2 COLUNAS E O GRÁFICO NA POSIÇÃO 4
imshow(rec_im);                                                            %EXIBE A IMAGEM RECUPERADA NA RECEPÇÃO DESTE PROCESSO                                     
title(sprintf('\\bfIMAGEM RECUPERADA\n\\rmBER: %.2g',ber));                %TÍTULO DO GRÁFICO PLOTADO

%CONFIGURAÇÕES EXTRAS RELACIONADAS A IMAGEM
set(gcf,'Position',[680 287 698 691]);                                     %CONFIGURAÇÃO DO POSICIONAMENTO DA FIGURA
save_file=0;                                                               %OPÇÃO PARA SALVAR O ARQUIVO

sPlotFig = scatterplot(SMHp2sfftsinal_final,1,0,'g.');                     %PLOTAGEM DO SINAL RECEBIDO COM RUÍDO
hold on                                                                    %RETENDO A PLOTAGEM FEITA
scatterplot(ytx,1,0,'k*',sPlotFig)                                         %PONTOS MQAM PLOTADOS NO GRÁFICO EM SOBREPOSIÇÃO AO DADO RECEBIDO

     end
         