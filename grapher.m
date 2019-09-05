testfolder="testsalpha0.2"
files=dir(testfolder)
for a=1:length(files)
    filename=files(a).name
    if filename~="." && filename~=".."  && filename~="figs"  
        filepath=fullfile(testfolder,filename)
        load(filepath)
        filename=char(filename)
        filename = filename(1:end-4)
        fig = figure
        plot(K,output)
        xlabel('K')
        ylabel('Proportion Found')
        title(filename)
        filepath = strcat(testfolder, "/figs/",filename,".fig")
        savefig(fig, filepath)
        close
    end
end