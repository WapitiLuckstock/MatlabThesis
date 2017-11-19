clear all
close all
clc
%% 
% excelfiles inladen
VerbruikInHalfUur=xlsread('VerbruiksprofielDH2003','sheet1','B2:AW366');  
TemperaturenInHalfUur=xlsread('Temperaturen2003','blad1','E2:AZ366');
RadiatieInHalfUur=xlsread('2016Instralingsflux.xls','blad1','B2:AW366');
TempdijleInKwartier=xlsread('tempdijle','Boortmeerbeek Weesbeek_','B9:B35048');

%gegevens en aannames
Verwarmingskost=20572.99; % doorgegeven via mail
Aantalapp=36; % niet zeker
Prijsgas=0.03280; % exclusief btw https://eneco.be/nl/elektriciteit-en-gas/prijs-gas
Verliezenoudnet=0.4; % aanname
Correctievoordrukkemomentenaanvoer=5; % +5-5 op aanvoer
Minbuiten=0;Maxbuiten=16;Minaanvoer=55;Maxaanvoer=70; % buitentemp -> aanvoertemp
warmtepompinout=18;
waterradgemoutput=20; %kwh
waterradplusmin=0.40;
exportgas=0.035; % aanname
prijswaterrad=5.500*waterradgemoutput; % aanname
totalekost=200000; %aanname
Areazonnecollector=10; % variabele aanname
ondergrensdijle = 3; % onder deze waarde neemt de warmtepomp geen dijle water op

% gegevenszonnecollector
%-> https://cdm.unfccc.int/Panels/ssc_wg/SSCWG7_repan3_Convers_fact_solar_collectors.pdf
%glazed flat plate
m0g=0.78;a1g=3.2;a2g=0.015;
%Evacuated Tube
m0e=0.76;a1e=1.2;a2e=0.008;
% gegevenspompenruben (aanname)
% afvoertemperatuur hangt tussen de 6 en de 14°C
Aanvoertemp=[45 55 70 75];
COP=[5.2 4.06 2.85 2.32];

% verband invoertemp dijle tov uitvoertemp met de COP ( gegevens van
% warmtepomp)
Aanvoertempdijle=[13.2 12.6 12 11.5 10.3 9.8 9.2 8.6];
Leavingtempdijle=7.5;
Deltatempdijle=Aanvoertempdijle-ones(size(Aanvoertempdijle))*Leavingtempdijle;
TheoretischeCOPwaarde=5.2;
COPdeltadijle=[5.2 5.27 5.4 5.63 5.12 4.34 3.41 2.55]/TheoretischeCOPwaarde;
b2=polyfit(Deltatempdijle,COPdeltadijle,2);
figure
grid on
X2=0:0.01:5.70;
Y2=b2(1)*X2.^2+b2(2)*X2+b2(3);
plot(X2,Y2)



%berekening van aangenomen verlies
[m,n]=size(VerbruikInHalfUur);
Aantalvakjes=m*n;
Voorbeeldverbruik=sum(sum(VerbruikInHalfUur));
Vraagperapp=(Verwarmingskost/Aantalapp)/Prijsgas;
Correctietoevoegingwarmte=Verliezenoudnet*Voorbeeldverbruik/Aantalvakjes;
Pereenheidcorrectie=Vraagperapp/(Voorbeeldverbruik+Correctietoevoegingwarmte*Aantalvakjes);
% omzetting van uurlijk totaalverbruik naar werkelijk uurlijk verbruik zonder verlies
Echtverbruik=VerbruikInHalfUur*Pereenheidcorrectie;

%correctiefactor op warmtevraag berekening
Gemiddeldewaardeperhalfuur=mean(Echtverbruik);
Maxgemtwaardeperhalfuur=max(Gemiddeldewaardeperhalfuur);
Mingemtwaardeperhalfuur=min(Gemiddeldewaardeperhalfuur);
[m2,n2]=size(Gemiddeldewaardeperhalfuur);
for i = 1:n2
    Gemiddeldewaardeperhalfuur(1,i)=(Gemiddeldewaardeperhalfuur(1,i)-Mingemtwaardeperhalfuur)/(Maxgemtwaardeperhalfuur-Mingemtwaardeperhalfuur)*10-5;
end

% buitentemp naar aanvoertemp
Aanvoertempvbuiten=TemperaturenInHalfUur;
[m3,n3]=size(Aanvoertempvbuiten);
for i = 1:m3
    for j = 1:n3
        Aanvoertempvbuiten(i,j)=Gemiddeldewaardeperhalfuur(1,j)+Maxaanvoer+(TemperaturenInHalfUur(i,j)*((Minaanvoer-Maxaanvoer)/(Maxbuiten-Minbuiten)));
        if Aanvoertempvbuiten(i,j)<55
            Aanvoertempvbuiten(i,j)=55;
        end
        if Aanvoertempvbuiten(i,j)>70
                Aanvoertempvbuiten(i,j)=70;
            
        end
    end
end

%% plot van de stooklijn
f1 = figure('Name','Stooklijn');
hplot1 = plot(Aanvoertemp,COP);
ax1 = gca;
grid on
hold on
b=polyfit(Aanvoertemp,COP,1); %y=ax+b -> b(1) = a, b(2) = b
X=1:80;
plot(X,b(1)*X+b(2))
xlabel('aanvoertemperatuur (°C)')
ylabel('COP')
title('Stooklijn')

%% Maken tabel watertemperaturen
tabeltempdijle=zeros(m,n);
for i = 1:m
    for j= 1:n
        tabeltempdijle(i,j)=(TempdijleInKwartier((((i-1)*2)*n)+(j*2))+TempdijleInKwartier((((i-1)*2)*n)+((j*2)-1)))/2;
        if tabeltempdijle(i,j) == 0
            tabeltempdijle(i,j)=tabeltempdijle(i-1,j);
        end
    end
end
minimumdijletemp=min(min(tabeltempdijle));
maximumdijletemp=max(max(tabeltempdijle));

for i = 1:m
    for j =1:n
        
    end
end

%% COP in functie van aanvoertemperatuur
[m4,n4]=size(Aanvoertempvbuiten);
COPwaarden=zeros(m,n);
factordijle=(b2(1)*(tabeltempdijle(i,j)-ondergrensdijle)^2+b2(2)*(tabeltempdijle(i,j)-ondergrensdijle)+b2(3));
for i = 1:m4
    for j= 1:n4
        if tabeltempdijle(i,j) < ondergrensdijle
            COPwaarden(i,j)=0;
        elseif tabeltempdijle(i,j) > max(Aanvoertempdijle)
            COPwaarden(i,j)=(b(1)*Aanvoertempvbuiten(i,j)+b(2));
        else
            COPwaarden(i,j)=(b(1)*Aanvoertempvbuiten(i,j)+b(2))*factordijle;
        end
    end
end
%% radiatie op warmtevraag
%warmtevoorziening zon : P/A=G*Mu-a1(Tm-Ta)-a2(Tm-Ta)² 
% Voor glazed
[m6,n6]=size(RadiatieInHalfUur);
KWhRadiatieg=zeros(m6,n6);
for i = 1:m6
    for j = 1:n6
        KWhRadiatieg(i,j)=Areazonnecollector*(RadiatieInHalfUur(i,j)*m0g-a1g*(Aanvoertempvbuiten(i,j)-TemperaturenInHalfUur(i,j))-a2g*(Aanvoertempvbuiten(i,j)-TemperaturenInHalfUur(i,j))^(2));
        if KWhRadiatieg(i,j)<0
            KWhRadiatieg(i,j)=0;
        end
    end
end

%radiatie effect op warmtevraag
% for i = 1:m6
%     for j = 1:n6
%         Echtverbruik(i,j)=Echtverbruik(i,j)-KWhRadiatieg(i,j);
%         if Echtverbruik(i,j)<0
%             Echtverbruik(i,j)=0;
%         end
%     end
% end

%% plotten warmtevraag duration curve
dcechtverbruik=zeros(365,24);
for i = 1:365
    for j = 1:24
        dcechtverbruik(i,j)=Echtverbruik(i,2*j)+Echtverbruik(i,2*j-1);
    end
end

maxEV=max(max(dcechtverbruik));
minEV=min(min(dcechtverbruik));
Aantalvaluesingrafiekwarmtevraag=round((maxEV-minEV)*100,0)+1;
q=0; % variabele voor volgende for lus
Grafiekwarmtevraag=zeros(Aantalvaluesingrafiekwarmtevraag,2);
for i = 1:Aantalvaluesingrafiekwarmtevraag
    for j = 1:365
        for k = 1:24
            if dcechtverbruik(j,k)>= q && dcechtverbruik(j,k) < (minEV+i/100)
                Grafiekwarmtevraag(i,2)=Grafiekwarmtevraag(i,2)+1;
            end
        end
    end
    q=minEV+i/100;
    Grafiekwarmtevraag(i,1)=(minEV+(i-1)/100)*Aantalapp; %warmtevraag voor het totale gebouw
end
[m5,~]=size(Grafiekwarmtevraag);

for i=1:(m5-1)
    Grafiekwarmtevraag(m5-i,2)=Grafiekwarmtevraag(m5-i,2)+Grafiekwarmtevraag((m5-(i-1)),2);
end
f2 = figure('Name','Warmtevraag');
hplot2= plot(Grafiekwarmtevraag(:,2),Grafiekwarmtevraag(:,1));
ax2 = gca;
grid on
xlabel('Tijd (uur)')
ylabel('Warmtevraag (kWh)')
title('Warmtevraag')
%% maken van elektriciteitsvraag voor de Duration Curve
Elektriciteitsvraag=zeros(m,n);
for i = 1:m
    for j = 1:n
        if COPwaarden(i,j)~=0
            Elektriciteitsvraag(i,j)=Echtverbruik(i,j)/COPwaarden(i,j);
        end
    end
end
%duration curve elektriciteitsvraag
dcelekvraag=zeros(365,24);
for i = 1:365
    for j = 1:24
        dcelekvraag(i,j)=Elektriciteitsvraag(i,2*j)+Elektriciteitsvraag(i,2*j-1);
    end
end
Maxelek=max(max(dcelekvraag));
Minelek=min(min(dcelekvraag));
Aantalvaluesingrafiekelektriciteitsvraag=round((Maxelek-Minelek)*100,0)+1;
q=0; % variabele voor volgende for lus
Grafiekelektriciteitsvraag=zeros(Aantalvaluesingrafiekelektriciteitsvraag,2);
for i = 1:Aantalvaluesingrafiekelektriciteitsvraag
    for j = 1:365
        for k = 1:24
            if dcelekvraag(j,k)>= q && dcelekvraag(j,k) < (Minelek+i/100)
                Grafiekelektriciteitsvraag(i,2)=Grafiekelektriciteitsvraag(i,2)+1;
            end
        end
    end
    q=Minelek+i/100;
    Grafiekelektriciteitsvraag(i,1)=(Minelek+(i-1)/100)*Aantalapp; %warmtevraag voor het totale gebouw
end
[m5,n5]=size(Grafiekelektriciteitsvraag);
GE=Grafiekelektriciteitsvraag; % zelfde tabel voor economisch model

for i=1:(m5-1)
    Grafiekelektriciteitsvraag(m5-i,2)=Grafiekelektriciteitsvraag(m5-i,2)+Grafiekelektriciteitsvraag((m5-(i-1)),2);
end
f3 = figure('Name','Elektriciteitsvraag');
plot(Grafiekelektriciteitsvraag(:,2),Grafiekelektriciteitsvraag(:,1))
ax3 = gca;
grid on
xlabel('Tijd (uur)')
ylabel('Elektriciteitsvraag (kWh)')

% Economisch model zonder zonneboiler en buffervat
[m5,~]=size(Grafiekelektriciteitsvraag);
Overschotelektriciteit=0;
Overelektabel=zeros(m5,2);
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && waterradgemoutput>=warmtepompinout
        Overelektabel(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabel(i,1)=Overelektabel(i,1)+waterradgemoutput-warmtepompinout;
        Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && waterradgemoutput>=warmtepompinout
        Overelektabel(i,1)=Overelektabel(i,1)+waterradgemoutput-warmtepompinout;
        Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
    elseif waterradgemoutput<warmtepompinout
        if waterradgemoutput>=Grafiekelektriciteitsvraag(i,1)
            Overelektabel(i,1)=waterradgemoutput-Grafiekelektriciteitsvraag(i,1);
            Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2);
            Overschotelektriciteit=Overschotelektriciteit+Overelektabel(i,1)*GE(i,2);
        else 
            Overelektabel(i,1)=0;
            Overelektabel(i,2)=Grafiekelektriciteitsvraag(i,2); 
        end
        
    end
end

% plus min afhankelijk van waterrad
Overelektabelmin=Overelektabel;
Overschotelektriciteitmin=0;
watermin=waterradgemoutput-waterradgemoutput*waterradplusmin;
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && watermin>=warmtepompinout
        Overelektabelmin(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmin(i,1)=Overelektabelmin(i,1)+watermin-warmtepompinout;
        Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && watermin>=warmtepompinout
        Overelektabelmin(i,1)=0;
        Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmin(i,1)=Overelektabelmin(i,1)+watermin-warmtepompinout;
        Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
    elseif watermin<warmtepompinout
        if watermin>=Grafiekelektriciteitsvraag(i,1)
            Overelektabelmin(i,1)=watermin-Grafiekelektriciteitsvraag(i,1);
            Overschotelektriciteitmin=Overschotelektriciteitmin+Overelektabelmin(i,1)*GE(i,2);
        else 
            Overelektabelmin(i,1)=0;
            Overelektabelmin(i,2)=Grafiekelektriciteitsvraag(i,2);
        end
    end
end
Overelektabelmax=Overelektabel;
Overschotelektriciteitmax=0;
watermax=waterradgemoutput+waterradgemoutput*waterradplusmin;
for i=1:m5
    if Grafiekelektriciteitsvraag(i,1)<=warmtepompinout && watermax>=warmtepompinout
        Overelektabelmax(i,1)=warmtepompinout-Grafiekelektriciteitsvraag(i,1);
        Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        Overelektabelmax(i,1)=Overelektabelmax(i,1)+(watermax-warmtepompinout);
        Overschotelektriciteitmax=Overschotelektriciteitmax+Overelektabelmax(i,1)*GE(i,2);
    elseif Grafiekelektriciteitsvraag(i,1)>warmtepompinout && watermax>=warmtepompinout
        Overelektabelmax(i,1)=0;
        Overelektabelmax(i,1)=Overelektabelmax(i,1)+(watermax-warmtepompinout);
        Overschotelektriciteitmax=Overschotelektriciteitmax+Overelektabelmax(i,1)*GE(i,2);
        
    elseif watermax<warmtepompinout
        if watermax>=Grafiekelektriciteitsvraag(i,1)
            Overelektabelmax(i,1)=watermax-Grafiekelektriciteitsvraag(i,1);
            Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        else 
            Overelektabelmax(i,1)=0;
            Overelektabelmax(i,2)=Grafiekelektriciteitsvraag(i,2);
        end
    end
end   
f4 = figure('Name','Elektriciteitsvraag en productie');
hplot4 = plot(Grafiekelektriciteitsvraag(:,2),Grafiekelektriciteitsvraag(:,1));
ax4 = gca;
grid on
hold on
plot(Overelektabel(:,2),Overelektabel(:,1))
hold on
plot(Overelektabelmin(:,2),Overelektabelmin(:,1))
hold on
plot(Overelektabelmax(:,2),Overelektabelmax(:,1))
xlabel(' Tijd (uren)')
ylabel(' Elektriciteitsvraag/overschot (kWh)')
title('Elektriciteitsoverschot voor export afhankelijk van productie waterrad')

f5 = figure('Name','Inkomsten');
kwhovertomoney=[ watermin Overschotelektriciteitmin*exportgas;waterradgemoutput Overschotelektriciteit*exportgas;watermax Overschotelektriciteitmax*exportgas];
hplot5 = plot(kwhovertomoney(:,1),kwhovertomoney(:,2));
ax5 = gca;
axis([10 30 0 10000])
grid on
xlabel('gemiddelde productie waterrad (kWh)')
ylabel('export inkomsten (euro)')
title('export inkomsten')


%% variabele kosten
extragas=0;
for i = 1:m
    for j = 1:n
        if Elektriciteitsvraag(i,j)*Aantalapp>warmtepompinout
            extragas=extragas+(Elektriciteitsvraag(i,j)*Aantalapp-warmtepompinout)*COPwaarden(i,j)*Prijsgas;
        end
    end
end
variabelekostperjaar= Overschotelektriciteit*exportgas+Verwarmingskost-extragas;
f6 = figure('Name','Variabele kosten');
Jaar=1:20;
hplot6 = plot(Jaar,variabelekostperjaar*Jaar);
ax6 = gca;
hold on
plot(Jaar,ones(size(Jaar))*totalekost)
grid on
xlabel('jaren')
ylabel('prijs')
title('terugverdientijd')

%% Waar zitten de pieken
Piekelek=Elektriciteitsvraag*Aantalapp;
for i = 1:m
    for j = 1:n
        if Elektriciteitsvraag(i,j)*Aantalapp<=warmtepompinout
            Piekelek(i,j)=0;
        end
    end
end
%hoe lang duren die pieken
Piekduur=zeros(m,3);
kollominpiekduur=1;
for i = 1:m-1
    for j=1:n-1
        if Piekelek(i,j+1)>0 && Piekelek(i,j)>0
          Piekduur(i,kollominpiekduur)=Piekduur(i,kollominpiekduur)+0.5;
        elseif Piekelek(i,j+1)==0 && Piekelek(i,j)>0
            kollominpiekduur=kollominpiekduur+1;
        end
    end
    kollominpiekduur=1;
end

% duration curve van de pieken
maxpiekduur=max(max(Piekduur));
minpiekduur=min(min(Piekduur));
durationpiek=zeros(maxpiekduur-minpiekduur,2);
for deltapiek=1:(maxpiekduur-minpiekduur)
    for i= 1:m
        for j=1:3
            if Piekduur(i,j)== deltapiek+minpiekduur
                durationpiek(deltapiek,2)=durationpiek(deltapiek,2)+0.5;
            end
            durationpiek(deltapiek,1)=minpiekduur+deltapiek;
        end
    end
end

for i=1:(maxpiekduur-minpiekduur-1)
    durationpiek(maxpiekduur-minpiekduur-i,2)=durationpiek(maxpiekduur-minpiekduur-i,2)+durationpiek((maxpiekduur-minpiekduur-(i-1)),2);
end
%% Maken van de Piekduur grafiek

f7 = figure('Name','Tekortkomingen');
hplot7 = plot(durationpiek(:,2),durationpiek(:,1));
ax7 = gca;
grid on
xlabel('aantal keren per jaar dat warmtepomp onvoldoende is')
ylabel('lengte van de piek (uur) ')
title('piekduur')
axis([0 60 0 8])

%% Plots samenvoegen op een overzichtelijk frame

f1.Visible = 'Off';
f2.Visible = 'Off';
f3.Visible = 'Off';
f4.Visible = 'Off';
f5.Visible = 'Off';
f6.Visible = 'Off';
f7.Visible = 'Off';
f8 = figure('Name','Totaaloverzicht','units','normalized','outerposition',[0 0 1 1]);
s1=subplot(4,2,1);
s2=subplot(4,2,2);
s3=subplot(4,2,3);
s4=subplot(4,2,4);
s5=subplot(4,2,5);
s6=subplot(4,2,6);
s7=subplot(4,2,7);
fig1=get(ax1,'children');
fig2=get(ax2,'children');
fig3=get(ax3,'children');
fig4=get(ax4,'children');
fig5=get(ax5,'children');
fig6=get(ax6,'children');
fig7=get(ax7,'children');
copyobj(fig1,s1);
copyobj(fig2,s2);  
copyobj(fig3,s3); 
copyobj(fig4,s4); 
copyobj(fig5,s5); 
copyobj(fig6,s6); 
copyobj(fig7,s7); 
disp('het einde');
%% te leveren warmte door gas

figure
% warmtevraag niet opgevangen door het waterrad
dcwnietopgevangen=dcechtverbruik*Aantalapp;
for i = 1:m
    for j = 1:n/2
        if (COPwaarden(i,2*j)+COPwaarden(i,2*j-1))/2~=0
            if dcelekvraag(i,j)*Aantalapp < warmtepompinout
                dcwnietopgevangen(i,j)=0;
            else
                dcwnietopgevangen(i,j)=dcwnietopgevangen(i,j)-warmtepompinout*(COPwaarden(i,2*j)+COPwaarden(i,2*j-1))/2;
            end
            if dcwnietopgevangen(i,j) <0
                dcwnietopgevangen(i,j)=0;
            end
        end
    end
end

Maxno=max(max(dcwnietopgevangen));
Minno=min(min(dcwnietopgevangen));
Aantalvaluesingrafiekwcnietopgevangen=round((Maxno-Minno),0)+1;
q=0; % variabele voor volgende for lus
Grafiekwnietopgevangen=zeros(Aantalvaluesingrafiekwcnietopgevangen,2);
for i = 1:Aantalvaluesingrafiekwcnietopgevangen
    for j = 1:365
        for k = 1:24
            if dcwnietopgevangen(j,k)>= q && dcwnietopgevangen(j,k) < (Minno+i)
                Grafiekwnietopgevangen(i,2)=Grafiekwnietopgevangen(i,2)+1;
            end
        end
    end
    q=Minno+i;
    Grafiekwnietopgevangen(i,1)=(Minno+(i-1)); %warmtevraag voor het totale gebouw
end

[m6,~]=size(Grafiekwnietopgevangen)
for i=1:(m6-1)
    Grafiekwnietopgevangen(m6-i,2)=Grafiekwnietopgevangen(m6-i,2)+Grafiekwnietopgevangen((m6-(i-1)),2);
end
plot(Grafiekwnietopgevangen(:,2),Grafiekwnietopgevangen(:,1))
grid on
axis([0 8760 0 max(Grafiekwnietopgevangen(:,1))])
hold on
plot(Grafiekwarmtevraag(:,2),Grafiekwarmtevraag(:,1))
Tebetalengasfactuur=0;
for i = 1:m
    for j =1:n/2
        Tebetalengasfactuur=Tebetalengasfactuur+dcwnietopgevangen(i,j);
    end
end
Tebetalengasfactuur=Tebetalengasfactuur*Prijsgas;





