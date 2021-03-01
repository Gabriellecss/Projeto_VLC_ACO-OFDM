% % Universidade do Estado do Rio de Janeiro - UERJ
% % Faculdade de Engenharia- FEN
% % Engenharia El�trica com �nfase em Telecomunica��es 
%% PROJETO DE CONCLUS�O DE CURSO 
%% Comunica��o por luz vis�vel (VLC):An�lise e s�ntese por simula��o utilizando o MATLAB
% % Gabrielle Cristina de Souza Silva
% % Leonardo Alexandre da Silva
%================================LOG================================
%�ltima atualiza��o - 01/02/2021
%==============================Descri��o============================
%Rotina para obter as curvas de BER x Eb/No do ACO-OFDM com Imagem
%===============================IN�CIO==============================

clear all;                                                                 %LIMPA AS VARIAVEIS DO CODIGO
close all;                                                                 %FECHA TODAS AS FIGURAS
format long;                                                               %AUMENTA O N�MERO DE CASAS DECIMAIS PARA MAIS DE 15 DIGITOS

% ACO-OFDM

%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%%                                                       +++++   PARAMETROS DE ENTRADA    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
z=1;                                                                       %DEFININDO VARI�VEL Z IGUAL A 4
while (z<10)                                                               %EXECUTA ENQUANTO A VARI�VEL Z FOR INFERIOR A 10
z=z+1;                                                                     %A VARI�VEL Z=Z+1
F=2^z;                                                                     %A VARI�VEL F RECEBE 2 ELEVADO AO VALOR DA VARI�VEL Z, NA PRIMEIRA RODAGEM=32

M = F;                                                                     %M � A ORDEM DA MODULA��O
k = log2(M);                                                               %K � O N�MERO DE BITS POR SIMBOLO
EbNoVec = (0:1:50)';                                                       %Eb/No VALORES QUE VARIA DE 0 A 15 COM ESPA�O DE 1 EM 1 (dB), REPRESENTADO EM 15 LINHAS
numSymPerFrame = 3000;                                                     %N�MERO DE QAM SIMBOLOS POR FRAME | !!COMO INTERFERE NA LINHA DO RESHAPE,TEM QUE SER DIVISIVEL POR 8
berEst = zeros(size(EbNoVec));                                             %MATRIX VAZIA PARA O BER ESTIMADO COM 15 LINHAS
cp=16;                                                                     %PREFIXO CICLICO

%_________________________________________________________________________________________________________________________
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
%% VARIA��O DE BER COM O SNR (dB)


for n = 1:length(EbNoVec)
    snrdB(n) = EbNoVec(n) + 10*log10(k);                                   %CONVERTENDO Eb/No to SNR em dB
    numErrs = 0;                                                           %DEFININDO A VARIAVEL DO NUMERO DE ERROS IGUAL A "0"
    numBits = 0;                                                           %DEFININDO A VARIAVEL DO N�MERO DE BITS IGUAL A "0"
    RR=0;
    
   %while numErrs < 600 
   while RR < 400 
     %******CONSIDERANDO UM DADO ALEAT�RIO COMO A ENTRADA DO SISTEMA*******
     RR=RR+1;                        
     data_source=randi([0 1],numSymPerFrame,k);                                 %GERADOR ALEAT�RIO EM BIN�RIO, COM NUMSYMPERFRAM LINHAS E K COLUNAS
     dataSym = bi2de(data_source);                                              %CONVERS�O DE BINARIO PARA DECIMAL, FORMANDO OS S�MBOLOS SORTEADOS

%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 1.                                                         +++++ GERANDO QAM    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
     qam_modulated_data =qammod(dataSym,M);                                     %GERANDO QAM (SIMBOLOS,ORDEM), SA�DA � S�O N�MEROS COMPLEXOS
     ytx =qam_modulated_data;                                                   %SA�DA DO QAM ytx=qam_modulated_data
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 2-3.                                                      +++++   SIMETRIA HERMITIANA    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
     
     %SERIAL PARA PARALELO
     ytxx=ytx.';                                                                 
     
     %HERMITIANA
     y = upsample(ytxx,2);                                                      %INSER��O DE ZEROS NAS SUBPORTADORAS PARES DO ACO-OFDM
     hermitian = [0 y conj(fliplr(y(1:length(y)-1)))];                          %APLICANDO HERMITIANO, ESPELHANDO A MATRIZ Y
     
     %IFFT
     ifftsinal=ifft(hermitian);                                                 %APLICANDO A IFFT NO SINAL, SAINDO N�MEROS REAIS EM COLUNAS NA MATRIZ
     
     %PARALELO PARA SERIAL
     p2sifftsinal=ifftsinal.';
     
     %ADICIONANDO PREFIXO CICLICO
     cpp2sifftsinal=[p2sifftsinal(end-cp+1:end,:);p2sifftsinal];

%   -----------------------------------------------------------------------------------------------------------------------------------------------------------
%% 8.                                           +++++   REALIZANDO O CORTE DA ONDA (CLIPPING)    +++++
%   -----------------------------------------------------------------------------------------------------------------------------------------------------------

%      %CLIPPING DO SINAL (CORTANDO VALORES ABAIXO DE ZERO)
     avg=0;                                                                %DEFININDO VARIAVEL AVG=0                                                                       
       clipped=cpp2sifftsinal;                                             %CLIPPED RECEBE OS VALORES DE 'OFDM_SINAL'
       for i=1:length(clipped)                                             %C�DIGO DE REPETI��O QUE VAI DE 1 AT� O COMPRIMENTO DE CLIPPED
 	     if clipped(i) > avg                                               %SE CLIPPED � MAIOR QUE AVG ENT�O
 		 clipped(i) = clipped(i);                                          %CLIPPED RECEBE O VALOR DO PROPRIO CLIPPED
     
         elseif clipped(i) < -avg                                          %SE CLIPPED < 0 ENT�O
 		    clipped(i) = 0;                                                %CLIPPED RECEBE O VALOR "0"
         end                                                               %ENCERRA O COMANDO END
       end                                                                 %ENCERRA O COMANDO FOR






       
     %SINAL RECEBIDO+RU�DO AWGN   
     sinalrecebido = awgn(clipped,snrdB(n),'measured');                    %SINAL RECEBIDO SENDO APLICADO O RU�DO AWGN
       
     %REMOVENDO PREFIXO CICLICO
     rmcpsinalrecebido=sinalrecebido(cp+1:end,:);                          %SINAL RECEBIDO SEM O PREFIXO CICLICO

     
     %SERIAL PARA PARALELO

     s2psinalrecebido=rmcpsinalrecebido.';                                 %MODIFICANDO O DADO PARA PARALELO
     
     %FFT
     
     fftsinal=fft(s2psinalrecebido);                                       %SINAL COM FFT
     
     %PARALELO PARA SERIAL
     
     p2sfftsinal=fftsinal.';                                               %DADO DE PARALELO PARA SERIAL
     
     %SIMETRIA HERMITIANA DEMOD

     SMHp2sfftsinal = p2sfftsinal(2:(length(p2sfftsinal)-1)/2+1);                             
     SMHp2sfftsinal_final = 2*downsample(SMHp2sfftsinal,2);
     
     
     %DEMODULA��O
     qamdemodulado = qamdemod(SMHp2sfftsinal_final,M);                     %DEMODULA O SINAL RECEBIDO, EM SERIAL, N�MEROS DECIMAIS
     dataOut = de2bi(qamdemodulado,'left-msb');                            %DADO RECEBIDO EM BIN�RIO
     dataout2=fliplr(dataOut);                                             %ESPELHANDO A MATRIZ RECEBIDA

     %ERROS DE BITS
     nErrors = biterr(data_source,dataout2);                               %CALCULO DO NUMERO DE ERRO DE BITS  
     
     numErrs = numErrs + nErrors;                                          %INCREMENTO DA QUANTIDADE DE ERRO DE BITS
     numBits = numBits + numSymPerFrame*k;                                 %INCREMENTO DA QUANTIDADE DE BITS ENVIADOS
     
   end                                                                     %FINALIZA��O DO COMANDO FOR
                                                                           %FINALIZA O LOOP DE ENQUANTO
                                                                           
     %C�LCULO DO BER ESTIMADO                                                                      
     berEst(n) = numErrs/numBits                                           %BER ESTIMADO = QUANTIDADE DE ERROS/QUANTIDADE DE BITS ENVIADOS
                                                                          
end
     %C�LCULO DO BER TE�RICO DO QAM COM RU�DO 
berTheory = berawgn(EbNoVec,'qam',M);                                      %C�LCULO  DO BER TE�RICO COM BASE NO Eb/N0 e a Modula��o QAM utilizada


%% *******GR�FICO FIGURA 7 - PLOTAGEM GR�FICO DE BER)*******

figure(1)                                                                  %FIGURA 7


title('ACO-OFDM Eb/No x BER 1024 subportadoras com 4,16,32,64,256,1024-QAM');


if M==4
semilogy(EbNoVec,berEst,'m-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')                   %LEGENDA DO GR�FICO
end


if M==16
semilogy(EbNoVec,berEst,'k-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')                   %LEGENDA DO GR�FICO
end


if M==32
semilogy(EbNoVec,berEst,'m-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')                   %LEGENDA DO GR�FICO
end

if M==64
semilogy(EbNoVec,berEst,'k-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')                   %LEGENDA DO GR�FICO
end

if M==128
semilogy(EbNoVec,berEst,'r-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')  %LEGENDA DO GR�FICO
end

if M==256
semilogy(EbNoVec,berEst,'b-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')  %LEGENDA DO GR�FICO
end

if M==1024
semilogy(EbNoVec,berEst,'r-o')                                             %PLOTAGEM DOS PONTOS * NO GR�FICO
hold on                                                                    %SEGURA O GR�FICO PARA SOBREPOR PLOTAGEM
grid on                                                                    %LINHAS DE GRADE DO GR�FICO
legend('4-QAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','1024-QAM')  %LEGENDA DO GR�FICO
end



%% C�LCULO DE BER TE�RICO E PLOTAGEM DO GR�FICO

%title('ACO-OFDM Eb/No x BER: Estimado vs Te�rico')
% semilogy(EbNoVec,berTheory,'*')                                          %PLOTAGEM DA CURVA TEORICA
% legend('BER 32-QAM','BER Te�rico')                                       %LEGENDA DO GR�FICO
% semilogy(EbNoVec,berEst,'r-*')                                           %PLOTAGEM DOS PONTOS * NO GR�FICO
% xlabel('Eb/No (dB)')                                                     %R�TULO DO EIXO X DO GR�FICO
% ylabel('Bit Error Rate M-QAM')                                           %R�TULO DO EIXO Y DO GR�FICO
% hold on

%ajustado 03-12-2020

end

