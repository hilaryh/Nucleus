function rescaledImg = fonf0(IMG)

for k = 1:size(IMG,1)
    IMG(k,:) = IMG(k,:)/IMG(k,1);
end

rescaledImg = IMG;