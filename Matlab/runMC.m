%%

seqs = cell(1,4);
seqs{1} = [repmat('A',[1 11]), repmat('C',[1 86])];
seqs{2} = [repmat('A',[1 12]), repmat('C',[1 85])];
seqs{3} = [repmat('A',[1 13]), repmat('C',[1 84])];
seqs{4} = [repmat('A',[1 14]), repmat('C',[1 83])];

pmcs = cell(1,4);

for i=1:4
    pmcs{i} = PoreMC(seqs{i},'samples',1e4);
end

for n=1:10000
    for i=1:4
        pmcs{i}.Next();
    end
    pmcs{4}.Plot();
    fprintf('\n');
end

%%

pmc = PoreMC(repmat('A',[1 27]),'samples',1e4,'V',.08);
for n=1:10000
    pmc.Next();
    if mod(n,10) == 0
        pmc.Plot();
    end
    fprintf('\n');
end