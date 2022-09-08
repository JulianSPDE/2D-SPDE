function createpeaks(vec)
    len = length(vec);
    for i=1:len
        SPDEsolvePeanutsinglepeak(100,400,1,vec(i),i)
    end
end
