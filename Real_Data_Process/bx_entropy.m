function E=bx_entropy(img)
%%%%%%%%%%%%%%%%ͼ����%%%%%%%%%%%
p=(abs(img)).^2;
PP=(sum(p));
S=PP/sum(PP);
E=-(sum(S.*log(S)));
end
    