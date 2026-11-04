function [order, deg] = vova_layer_order(H, Zc, name)
%VOVA_LAYER_ORDER Порядок обхода слоёв для VOVA_LAYERED_BP_DECODE.
%   [order, deg] = VOVA_LAYER_ORDER(H, Zc, name)
%     name: 'forward'     - прямой, 1 .. nLayers
%           'reverse'     - обратный, nLayers .. 1
%           'degree'      - по степени строк, от большей к меньшей
%                           (при равной степени - прямой порядок)
%           'degree-asc'  - по степени строк, от меньшей к большей
%   deg - степень строки каждого слоя (число единиц в строке H)

nLayers = size(H, 1) / Zc;
deg = zeros(1, nLayers);
for l = 1:nLayers
    deg(l) = nnz(H((l - 1) * Zc + 1, :));
end
switch lower(name)
    case 'forward',    order = 1:nLayers;
    case 'reverse',    order = nLayers:-1:1;
    case 'degree',     [~, order] = sort(deg, 'descend');
    case 'degree-asc', [~, order] = sort(deg, 'ascend');
    otherwise, error('vova_layer_order:name', 'Неизвестный порядок: %s', name);
end
order = order(:).';
end
