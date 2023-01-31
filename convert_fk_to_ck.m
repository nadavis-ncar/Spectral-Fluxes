function [K_e_new, K_w_new] = convert_fk_to_ck(freq, wavenum, cos_lat, c, K_e, K_w, time_res, lon_fraction)

RADIUS=6.371e6;          
DAY=86400/time_res;           
r_cos1=RADIUS*cos_lat*lon_fraction;       
nk=length(wavenum);      
nf=length(freq);     
nc=length(c);
df=(2*pi)/(2*DAY)/(nf-1); 

[k_axis1,f_axis1]=meshgrid(wavenum,freq);
[k_axis2,c_axis2]=meshgrid(wavenum,c);
f_axis2 = c_axis2.*(k_axis2/r_cos1)/(2*pi)*DAY;

K_e_new=griddata(k_axis1,f_axis1,squeeze(K_e)/df,k_axis2,double(f_axis2)).*(k_axis2/r_cos1);
K_w_new=griddata(k_axis1,f_axis1,squeeze(K_w)/df,k_axis2,double(f_axis2)).*(k_axis2/r_cos1);


for k=1:nk
for cc=1:nc
  dc=df*r_cos1/wavenum(k); % unresolved waves
  if(c(cc)<dc)
    K_e_new(cc,k)=nan;
    K_w_new(cc,k)=nan;
  else
    if(isnan(K_e_new(cc,k)))
      K_e_new(cc,k)=0;
    end
    if(isnan(K_w_new(cc,k)))
      K_w_new(cc,k)=0;
    end
  end  
end
end
return
