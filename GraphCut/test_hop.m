load('./tmp_hop');
imsize = size(imgs{1});
% test hop
nl = size(Dc,1);
Sc = ones(nl)-diag(ones(1,nl));
% gch = GraphCut('open_hop', Dc, Sc, sG, hop);
% [gch L_hop] = GraphCut('expand_QPBO',gch);
% L_hop = reshape(double(L_hop)+1,imsize(1:2));
% [gch se_hop de_hop he_hop] = GraphCut('energy', gch);
% % [gch se_hop de_hop] = GraphCut('energy', gch);
% gch = GraphCut('close', gch);

gch = GraphCut('open_hop', Dc, Sc, sG, hop_gen, hop);
[gch Lgc_hop] = GraphCut('expand',gch);
Lgc_hop = reshape(double(Lgc_hop)+1,imsize(1:2));
[gch segc_hop degc_hop hegc_hop] = GraphCut('energy', gch);
gch = GraphCut('close', gch);

% gch = GraphCut('open', Dc, Sc, sG);
% [gch L] = GraphCut('expand_QPBO',gch);
% L = reshape(double(L)+1,imsize(1:2));
% [gch se_hop de_hop] = GraphCut('energy', gch);
% gch = GraphCut('close', gch);

gch = GraphCut('open', Dc, Sc, sG);
[gch Lgc] = GraphCut('expand',gch);
Lgc = reshape(double(Lgc)+1,imsize(1:2));
[gch segc degc] = GraphCut('energy', gch);
gch = GraphCut('close', gch);

figure(1);
% subplot(2,2,1);
% imagesc(L_hop); title('L_hop');
subplot(2,2,2);
imagesc(Lgc_hop); title('Lgc_hop');
% subplot(2,2,3);
% imagesc(L); title('Lgc');
subplot(2,2,4);
imagesc(Lgc); title('Lgc');