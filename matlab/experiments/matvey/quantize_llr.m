function q = quantize_llr(llr, nBits, Lsat)
%QUANTIZE_LLR Равномерное симметричное квантование LLR с насыщением.
%   q = QUANTIZE_LLR(llr, nBits, Lsat)
%   Порт quantize_llr из task3_llr_quantization.m
%   (https://github.com/MobinKhatib/information-theory-ldpc-project):
%   LLR обрезается до [-Lsat, Lsat] и округляется до одного из 2^nBits - 1
%   уровней -Qmax..Qmax, Qmax = 2^(nBits-1) - 1, шаг Lsat/Qmax.
%   Возвращаются восстановленные вещественные значения (int * шаг), чтобы
%   их можно было подать в nrLDPCDecode. nBits = Inf - без квантования.

if isinf(nBits)
    q = llr;
    return
end
Qmax = 2^(nBits - 1) - 1;
step = Lsat / Qmax;
q = round(min(max(llr, -Lsat), Lsat) / step) * step;
end
