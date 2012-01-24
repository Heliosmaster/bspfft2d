clear all; close all;
%for n=[128 256 512 1024 2048]
    n2 = 2048;
    n1 = 1024;
    str1 = [int2str(n1)];
    str = [ str1 '.txt'];
    f = dlmread(str);
    trasp = f(:,[1 2]);
    mxn  = f(:,[1 3]);

    trasp_1 = trasp(1,2);
    mxn_1 = mxn(1,2);

    u = 1;
    for i=2:size(trasp,1)
        r = trasp_1/trasp(i,2);
        u = [u r];
    end
    
    v = 1;
    for i=2:size(mxn,1)
        s = mxn_1/mxn(i,2);
        v = [v s];
    end
    
    plot(trasp(1:6,1),trasp(1:6,1),'r');
    hold on;
    plot(trasp(1:6,1),u(1:6),'-*');
    plot(trasp(1:6,1),v(1:6),'--og');
    legend('Ideal','Row cyclic','M \times N');
    legend('Location','NorthWest');
    tit = ['Speedup plot for a ' str1 ' matrix'];
    title(tit);
    xlabel('p');
    ylabel('Speedup');
    set(gcf,'Position',[386 91 588 567]);
    print(gcf,'-depsc','-r300',str1);
    close all;
%end