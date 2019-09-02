files=dir("tests")
for a=1:length(files)
    filename=files(a).name
    if filename~="." && filename~=".."  && filename~="figs"  
        filepath=fullfile("tests",filename)
        load(filepath)
        filename=char(filename)
        filename = filename(1:end-4)
        fig = figure
        plot(K,output)
        xlabel('K')
        ylabel('Proportion Found')
        title(filename)
        filepath = strcat("tests/figs/",filename,".fig")
        savefig(fig, filepath)
        close
    end
end