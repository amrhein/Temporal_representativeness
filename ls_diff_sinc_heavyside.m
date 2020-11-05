function [out] = ls_diff_sinc_heavyside(t,rtp)

fun1 =  @(nu,tp) (1-sinc(tp*nu)).^2;
fun2 =  @(nu,tp) (sinc(tp*nu)).^2;

uplim = 1;
%rtp = 1:10:2000;

for ii = 1:length(rtp)
out(ii) = integral(@(nu)fun1(nu,rtp(ii)),0  ,1/(2*t)) + ...
          integral(@(nu)fun2(nu,rtp(ii)),1/(2*t),uplim);
end
%out = [];

out = sqrt(out);