function dt=str2datetime(str)
%function dt=str2datetime(str)
%
% Groppe hack

[deight, tyme]=strtok(str,' ');
deight_bkdwn=textscan(deight,'%d-%d-%d');
Y=deight_bkdwn{1};
M=deight_bkdwn{2};
D=deight_bkdwn{3};

tyme_bkdwn=textscan(tyme,' %d:%d:%f');
H=tyme_bkdwn{1};
MI=tyme_bkdwn{2};
S=tyme_bkdwn{3};

dt=datetime(Y,M,D,H,MI,S);
