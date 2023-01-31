% Space-Time Spectral Analysis
%
% Computing the space-time spectrum for the covariance of two time series
% Format: [East,West] = space_time(u,v)
%
% Dimensions of the input (u,v)
% 2d: longitude-time
% 3d: longitude-time-latitude
%
% Dimensions of the output (East, West)
% 2d: wavenumber-frequency
% 3d: wavenumber-frequency-latitude
%
% East (West) denotes the eastward (westward) propagating waves.

% Note: The time spectrum of a real time series is symmetric about the nyquist frequency, 
% so we only take the first half of the spectrum in the end.

function [East,West] = space_time_new(u,v,windowing)
warning off MATLAB:divideByZero

%For limited-longitude analysis
if strcmp(windowing,'true')
   ham=0.54-0.46*cos(2*pi*[0:size(u,1)-1]/(size(u,1)-1));
   u=u.*repmat(ham(:),[1 size(u,2) size(u,3)]);
   v=v.*repmat(ham(:),[1 size(u,2) size(u,3)]);
end

if(ndims(u) == 2)
   [East,West] = space_time_2d(u,v);
elseif(ndims(u) == 3)
   [East,West] = space_time_3d(u,v);
end

return

%=================================================================================
function [East,West] = space_time_2d(u_xt,v_xt)

num_lon=size(u_xt,1);
num_time=size(u_xt,2);

%-----

u_kt = fft(u_xt,[],1)/num_lon;
C_u_kt = real(u_kt);
S_u_kt = -imag(u_kt);

C_u_kw = fft(C_u_kt,[],2)/num_time;
A_u = real(C_u_kw);
B_u = -imag(C_u_kw);

S_u_kw = fft(S_u_kt,[],2)/num_time;
a_u = real(S_u_kw);
b_u = -imag(S_u_kw);

K_u_kw = (A_u-b_u).^2+(-B_u-a_u).^2; % k,w
K_u_k_w= (A_u+b_u).^2+(B_u-a_u).^2; % k,-w
phi_u_kw = atan((-B_u-a_u)./(A_u-b_u));
phi_u_k_w= atan((+B_u-a_u)./(A_u+b_u));

%------
v_kt = fft(v_xt,[],1)/num_lon;
C_v_kt = real(v_kt);
S_v_kt = -imag(v_kt);

C_v_kw = fft(C_v_kt,[],2)/num_time;
A_v = real(C_v_kw);
B_v = -imag(C_v_kw);

S_v_kw = fft(S_v_kt,[],2)/num_time;
a_v = real(S_v_kw);
b_v = -imag(S_v_kw);

K_v_kw = (A_v-b_v).^2+(-B_v-a_v).^2; % k,w
K_v_k_w= (A_v+b_v).^2+(B_v-a_v).^2; % k,-w
phi_v_kw = atan((-B_v-a_v)./(A_v-b_v));
phi_v_k_w= atan((+B_v-a_v)./(A_v+b_v));

%------
W = (((A_u-b_u).*(A_v-b_v)+(-B_u-a_u).*(-B_v-a_v)))*2;
E = (((A_u+b_u).*(A_v+b_v)+(+B_u-a_u).*(+B_v-a_v)))*2;
West(1:(num_lon/2+1),1:(num_time/2+1))=W(1:(num_lon/2+1),1:(num_time/2+1));
East(1:(num_lon/2+1),1:(num_time/2+1))=E(1:(num_lon/2+1),1:(num_time/2+1));

return

%=================================================================================
function [East,West] = space_time_3d(u_xty,v_xty)

num_lon=size(u_xty,1);
num_time=size(u_xty,2);

%-----
u_kty = fft(u_xty,[],1)/num_lon;
C_u_kty = real(u_kty);
S_u_kty = -imag(u_kty);

C_u_kwy = fft(C_u_kty,[],2)/num_time;
A_u = real(C_u_kwy);
B_u = -imag(C_u_kwy);

S_u_kwy = fft(S_u_kty,[],2)/num_time;
a_u = real(S_u_kwy);
b_u = -imag(S_u_kwy);

K_u_kwy = (A_u-b_u).^2+(-B_u-a_u).^2; % k,w
K_u_k_wy= (A_u+b_u).^2+(B_u-a_u).^2; % k,-w
phi_u_kwy = atan((-B_u-a_u)./(A_u-b_u));
phi_u_k_wy= atan((+B_u-a_u)./(A_u+b_u));

%------
v_kty = fft(v_xty,[],1)/num_lon;
C_v_kty = real(v_kty);
S_v_kty = -imag(v_kty);

C_v_kwy = fft(C_v_kty,[],2)/num_time;
A_v = real(C_v_kwy);
B_v = -imag(C_v_kwy);

S_v_kwy = fft(S_v_kty,[],2)/num_time;
a_v = real(S_v_kwy);
b_v = -imag(S_v_kwy);

K_v_kwy = (A_v-b_v).^2+(-B_v-a_v).^2; % k,w
K_v_k_wy= (A_v+b_v).^2+(B_v-a_v).^2; % k,-w
phi_v_kwy = atan((-B_v-a_v)./(A_v-b_v));
phi_v_k_wy= atan((+B_v-a_v)./(A_v+b_v));

%------
W = (((A_u-b_u).*(A_v-b_v)+(-B_u-a_u).*(-B_v-a_v)))*2;
E = (((A_u+b_u).*(A_v+b_v)+(+B_u-a_u).*(+B_v-a_v)))*2;
West(1:(num_lon/2+1),1:(num_time/2+1),:)=W(1:(num_lon/2+1),1:(num_time/2+1),:);
East(1:(num_lon/2+1),1:(num_time/2+1),:)=E(1:(num_lon/2+1),1:(num_time/2+1),:);

return
