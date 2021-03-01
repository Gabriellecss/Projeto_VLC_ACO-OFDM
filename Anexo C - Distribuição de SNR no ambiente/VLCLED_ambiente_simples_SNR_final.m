% % %Referencia do codigo de  Er. PANCHAL PRATIK
% % %https://www.mathworks.com/matlabcentral/answers/335713-simple-code-for-vlc
%%  Universidade do Estado do Rio de Janeiro - UERJ
% % Faculdade de Engenharia- FEN
% % Engenharia El�trica com �nfase em Telecomunica��es 
%% PROJETO DE CONCLUS�O DE CURSO 
%% Comunica��o por luz vis�vel (VLC):An�lise e s�ntese por simula��o utilizando o MATLAB
% % Gabrielle Cristina de Souza Silva
% % Leonardo Alexandre da Silva
%================================LOG================================
% �ltima atualiza��o - 21/02/2021
%==============================Descri��o============================
%% Rotina de performance do SNR em ambiente INDOOR - Para sistema VLC %%
%===============================IN�CIO==============================
clear all;
clc;
close all;

%% % PARAMETROS B�SICOS REQUERIDOS %
Incidence = 70*pi/180;
TX_FOV = 70;                % Campo de vis�o (FOV) do transmissor
RX_FOV = 90;                % Campo de vis�o (FOV) do receptor
Tx = [4,4,4];               % Localiza��o do transmissor - LED
% Rxp = [2,2];                % Localiza��o do receptor - Fotodiodo 
W_Room = 4;                 % Ambiente - Largura da Sala
L_Room = 4;                 % Ambiente - Comprimento da Sala
H_Room = 3;                 % Ambiente - Altura entre o transmissor e o receptor
R = 1;                      % Responsividade do fotodiodo
Apd = 1e-4;                 % �rea do fotodiodo
Rb = 1e6;                   % Taxa de dados do sistema
Iamp = 5e-12;               % Corrente do amplificador
q = 1.6e-19;                % Carga do eletron
Bn = 50e6;                  % Largura de banda de ru�do
I2 = 0.562;                 % Fator de largura de banda de ru�do
PLED = 1;                   % Pot�ncia emitida por LED (W)
index =1;
HLED = 1;
%Configura��o do plano e dos objetos
[W L] = meshgrid(-(W_Room/2) : 0.10 : (W_Room/2));      % Considere o comprimento do bloco para a sala
xydist = sqrt((W).^2 + (L).^2);
hdist = sqrt(xydist.^2 + HLED.^2);
%D = Tx - Rx;
%d = norm(D);
%Incidence = acos()
A_Irradiance = ((Tx(3)-HLED)./hdist);
%I(index) = Irradiance*180/pi;
%if abs(Incidence <= RX_FOV)
      p = TX_FOV ;
      Tx_FOV = (TX_FOV*pi)/180;
      %%  C�LCULO B�SICO NO SISTEMA VLC %
      % Parametro Lambertiano 
      m = real(-log(2)/log(cos(Tx_FOV))); 
      % Intensidade de radia��o em um ponto particular
      Ro = real(((m+1)/(2*pi)).*A_Irradiance^m);
      % Pot�ncia transmitida por LED
      Ptx = PLED .* Ro;
      %% Ganho de canal (coeficiente de canal do canal LOS)
      
      %Theta=atand(sqrt(sum((Tx-Rx).^2))/H_Room);
      HLOS = (Apd./hdist.^2).*cos(Incidence).*Ro;
      
      % Energia recebida por fotodiodo
      Prx = HLOS.*Ptx;
      % Calcule o ru�do no sistema
      Bs = Rb*I2;
      Pn = Iamp/Rb;
      Ptotal = Prx+Pn;
      new_shot = 2*q*Ptotal*Bs;
      new_amp = Iamp^2*Bn;
      % Calcular SNR
      new_total = new_shot + new_amp;
      %SNRl = 100;teste
      SNRl = (R.*Prx).^2./ new_total;
      SNRdb = 10*log10(SNRl);

%             else
%                     SNRl = 0;
%                     SNRdb = 0;
%             end
index = index + 1;

%% Gr�fico de plotagem %

%%Ambiente em 3D - SNR
figure(1);
surf(W,L,SNRdb);
%%Outros tipos de plot
% meshc(W,L,SNRdb) % 3D + 2D mesmo plot
% mesh(SNRdb); % 3D sem preenchimento
%ylim([0 30]);
title('Distribui��o 3D do SNR na sala');
xlabel('Comprimento da Sala(metros)');
ylabel('Largura da Sala (metros)');
zlabel('SNR em dB');
% Escala de sinal- colorbar
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'Nivel de Sinal (db)';
set(colorTitleHandle ,'String',titleString);

%Ambiente em 2D - Contorno
figure(2);
contourf(W,L,SNRdb);
title('Distribui��o 2D do SNR na sala - Contorno 2D');
xlabel('Comprimento da Sala (metros)');
ylabel('Largura da Sala(metros)');
zlabel('SNR em dB');
% Escala de sinal- colorbar
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'Nivel de Sinal (db)';
set(colorTitleHandle ,'String',titleString);
