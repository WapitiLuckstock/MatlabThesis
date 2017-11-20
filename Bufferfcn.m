function GebufferdeWarmtevraag = Bufferfcn(Warmtevraag,Warmtepompvermogen,COP)
c_p = 4.18; %kJ/kgK warmtecoefficient van water
V_tank = 2000; %Tank Volume in Liter
t_filled = V_tank/fillrate; %tijd tot de tank volledig gevuld is
t=0:1:t_filled;
water_volume_in_tank = fillrate*t;

T_storage_winter = 273+70;
T_discharged_winter = 273+50;

Warmtevraaginhalfuur = Warmtevraag;
Warmtevraaginuur = turnHalfHourtoHourBySumming(Warmtevraaginhalfuur);
Warmtepompvermogen = Warmtepompvermogen; %Productieprofiel van warmtepomp, deze gaat nog nauwkeuriger worden met waterradproductieprofiel
COPmat = COP; %berekende sCOP
W_water = water_volume_in_tank*c_p*(T_storage_winter-T_discharged_winter)*(1/3600); %opgeslagen warmte in kWh
W_tank = max(W_water); %maximaal opgslaagbare warmte in de tank

Warmtepompwarmte = ones(365,24);
for i = 1:24
    for j = 1:365
           Warmtepompwarmte(j,i)=Warmtepompvermogen(j,i)*COP(j,i);
    end
end
%Berekening van hoeveelheden kWh tekort er is om aan warmtevraag te voldoen
Warmtepomptekort = ones(365,24);
for i = 1:24
    for j = 1:365
   Warmtepomptekort(j,i)=Warmtevraaginuur(j,i)-Warmtepompwarmte(j,i);
   if(Warmtepomptekort(j,i)<0)
       Warmtepomptekort(j,i)=0;
   end
    end
end
%Berekening van hoeveelheden kWh er teveel zijn om aan de warmtevraag te
%voldoen, deze zijn opslaagbaar in de buffer
Warmtepompteveel = ones(365,24);
for i = 1:24
    for j = 1:365
        Warmtepompteveel(j,i) = Warmtepompwarmte(j,i)-Warmtevraaginuur(j,i);
        if(Warmtepompteveel(j,i)<0)
            Warmtepompteveel(j,i) = 0;
        end
    end
end
%Het maken van de buffermatrices
Bufferopslag = zeros(365,24);
Warmtevraagmetbuffer = zeros(365,24);
Bufferoverschot = zeros(365,24);
for j = 1:365
    for i = 1:24
        if(i==1 && j == 1) %Opstart van Bufferopslag
        Bufferopslag(j,i) = Warmtepompteveel(j,i);
        elseif (i>1) %Uren na opstart
        Bufferopslag(j,i) = Warmtepompteveel(j,i)+Bufferoverschot(j,i-1);
        elseif (i == 1 && j>1) %Dag verder
        Bufferopslag(j,i) = Bufferopslag(j-1,24);        
        end
        
        if(Bufferopslag(j,i)>W_tank) %Begrenzing van max op te slagen energie
        Bufferopslag(j,i)=W_tank;
        end  
        
        if(Warmtevraaginuur(j,i)>Warmtepompwarmte(j,i))%Ontladen bij hogere warmtevraag dan levering
        Warmtevraagmetbuffer(j,i) = Warmtevraaginuur(j,i)-Bufferopslag(j,i); %Volledige bufferdepletie
            if(Warmtevraagmetbuffer(j,i)<Warmtepompwarmte(j,i)) %Bufferdepletie als niet volledige buffer leeg geraakt
            Bufferoverschot(j,i) = Warmtepompwarmte(j,i)-Warmtevraagmetbuffer(j,i); %Berekening van overblijvende buffer
            Warmtevraagmetbuffer(j,i) = Warmtepompwarmte(j,i); %Gecorrigeerde warmtevraag met buffer
            end
        elseif(Warmtevraaginuur(j,i)<=Warmtepompwarmte(j,i)) %Opladen Buffer
            Warmtevraagmetbuffer(j,i) = Warmtevraaginuur(j,i); %Ogenblikkelijke warmtevraag
            Bufferoverschot(j,i) = Bufferopslag(j,i);
        end
    end
end
end