%histogram data collection

width_0b = 1000*[-0.0256,-0.023055,-0.02051,-0.017965,-0.01542,-0.012875,-0.01033,-0.007785,-0.00524,-0.002695,-0.00015,0.002395,0.00494,0.007485,0.01003,0.012575,0.01512,0.017665,0.02021,0.022755,0.0253];
bin_0b = [13,4,51,152,275,417,702,1167,1838,1658,1780,1772,1321,924,735,466,331,233,127,31];

width_10b = 1000*[-0.0275,-0.024815,-0.02213,-0.019445,-0.01676,-0.014075,-0.01139,-0.008705,-0.00602,-0.003335,-0.00065,0.002035,0.00472,0.007405,0.01009,0.012775,0.01546,0.018145,0.02083,0.023515,0.0262];
bin_10b = [10,11,74,137,283,362,624,1105,1633,1803,1511,1727,1172,917,588,290,259,117,76,1];

width_0i = 1000*[-0.0322,-0.028095,-0.02399,-0.019885,-0.01578,-0.011675,-0.00757,-0.003465,0.00064,0.004745,0.00885,0.012955,0.01706,0.021165,0.02527,0.029375,0.03348,0.037585,0.04169,0.045795,0.0499];
bin_0i = [35,72,96,127,276,613,1606,1985,1853,1525,966,714,690,506,341,181,115,51,71,35];

width_10i = 1000*[-0.0397,-0.033675,-0.02765,-0.021625,-0.0156,-0.009575,-0.00355,0.002475,0.0085,0.014525,0.02055,0.026575,0.0326,0.038625,0.04465,0.05067,0.0567,0.062725,0.06875,0.074775,0.0808];
bin_10i = [59,89,148,300,391,749,1236,1095,864,517,365,134,59,20,0,1,0,1,0,1];

width_10b_lo = 1000*[-0.0285,-0.02577,-0.02304,-0.02031,-0.01758,-0.01485,-0.01212,-0.00939,-0.00666,-0.00393,-0.0012,0.00153,0.00426,0.00699,0.00972,0.01245,0.01518,0.01791,0.02064,0.02337,0.0261];
bin_10b_lo = [1,14,24,64,145,238,344,753,1095,1541,1401,1498,1514,1164,727,447,268,125,79,31];

width_10i_lo = 1000*[-0.031,-0.026985,-0.02297,-0.018955,-0.01494,-0.010925,-0.00691,-0.002895,0.00112,0.005132,0.00915,0.013165,0.01718,0.021195,0.02521,0.029225,0.03324,0.037255,0.04127,0.045285,0.0493];
bin_10i_lo = [43,42,51,73,73,242,710,953,998,842,649,430,424,362,268,182,65,29,20,4];

width_3deg_b = 1000*[-0.0247,-0.02203,-0.01936,-0.01669,-0.01402,-0.01135,-0.00868,-0.00601,-0.00334,-0.00067,0.002,0.00467,0.00734,0.01001,0.01268,0.01535,0.01802,0.02069,0.02336,0.02603,0.0287];
bin_3deg_b = [22,35,89,293,530,807,1468,2273,2402,2048,2340,1747,1381,1134,770,434,333,95,57,38];

width_3deg_i = 1000*[-0.0323,-0.028445,-0.02459,-0.020735,-0.01688,-0.013025,-0.00917,-0.005315,-0.00146,0.002395,0.00625,0.010105,0.01396,0.017815,0.02167,0.025525,0.02938,0.033235,0.03709,0.040945,0.0448];
bin_3deg_i = [9,31,11,37,27,77,370,657,618,597,537,294,282,256,186,175,138,26,10,27];

width_10b_ir = 1000*[-0.0288, -0.02521,-0.02162,-0.01803,-0.01444,-0.01085,-0.00726,-0.00367,-0.00008, 0.00351,0.0071,0.01069,0.01428,0.01787,0.02146,0.02505,0.02864,0.03223,0.03582,0.03941,0.043];
bin_10b_ir = [14, 61, 197, 519, 778, 1680, 2812, 3429, 3306, 2826, 2276, 928, 496, 242, 59, 2, 0, 1,1,0];

width_10i_ir = 1000*[-0.0312, -0.027175, -0.02315, -0.019125, -0.0151, -0.011075, -0.00705, -0.003025, 0.001, 0.005025, 0.00905, 0.013075, 0.0171, 0.021125, 0.02515, 0.029175, 0.0332, 0.037225, 0.04125, 0.0452750,0.0493];
bin_10i_ir = [313, 312, 391, 468, 468, 762, 1441, 1536, 1324, 865, 557, 545, 491, 351, 369, 261, 45, 46, 25,19];

figure
hold on
subplot(5,2,1)
histogram('BinEdges',width_0b,'BinCount',bin_0b,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('0^{\circ} Simulated Bennu, \mu = 1.198 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,3)
histogram('BinEdges',width_10b,'BinCount',bin_10b,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('10^{\circ} Simulated Bennu, \mu = 0.167 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,2)
histogram('BinEdges',width_0i,'BinCount',bin_0i,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('0^{\circ} Simulated Itokawa, \mu = 4.941 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,4)
histogram('BinEdges',width_10i,'BinCount',bin_10i,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('10^{\circ} Simulated Itokawa, \mu = 2.573 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,5)
histogram('BinEdges',width_10b_lo,'BinCount',bin_10b_lo,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Bennu - Limb Only, \mu = 1.701 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,6)
histogram('BinEdges',width_10i_lo,'BinCount',bin_10i_lo,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Itokawa - Limb Only, \mu = 6.939 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,7)
histogram('BinEdges',width_3deg_b,'BinCount',bin_3deg_b,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('3^{\circ} Images - Bennu, \mu = 1.370 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,8)
histogram('BinEdges',width_3deg_i,'BinCount',bin_3deg_i,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('3^{\circ} Images - Itokawa, \mu = 6.417 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,9)
histogram('BinEdges',width_10b_ir,'BinCount',bin_10b_ir,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Bennu - Infrared Fusion, \mu = 0.318 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on
subplot(5,2,10)
histogram('BinEdges',width_10i_ir,'BinCount',bin_10i_ir,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
xlabel('Error (m)','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Itokawa - Infrared Fusion, \mu = 0.997 m','FontSize',14)
xlim([-40 40])
ylim([0 3500])
grid on


xlim([-40 40])
ylim([0 3500])
