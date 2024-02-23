clear;
clc;
cursor =[3914,12048,15972,24074,27952,36064,39986,48104,51997,60050,63956,72023,75964,84071,88004,96042];
block = zeros(8,2);
for i = 1:2:15
    block((i+1)/2,1) = cursor(i);
    block((i+1)/2,2) = cursor(i+1);
end

for k = 1:8
    [x,fs] = audioread("DTMF_dialing.ogg",[block(k,1),block(k,2)]);
    sound(x,fs);
    pause(1);
    mark = linspace(block(k,1),block(k,2),11);
    numbers = zeros(1,10);
    for j = 1:10
        [fr,fc] = dialing_frequency(mark(j),mark(j+1));
        numbers(j) = find_number(fr,fc);
    end
    display(['第',num2str(k),'段数字为：',num2str(numbers)])
    pause(0.5)
end


function [fr,fc] = dialing_frequency(t1,t2)
    % input: 输入音频的开始和结束时间（音频中只包含【一个】拨号音）
    % ouput: 返回低频fr和高频fc
    [x,fs] = audioread("DTMF_dialing.ogg",[round(t1),round(t2)]); % round是确保t1,t2为整数值
    n = length(x);
    X = fft(x);
    Y = fftshift(X);
    powershift = abs(Y).^2/n;     % zero-centered power
    fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
    fshift = fshift(:);
%     plot(fshift,powershift);
    zero_mark = round(length(fshift)/2);
    for i = zero_mark:length(fshift)
        if fshift(i)-1000 < 5
            one_thousand_mark = i;
        end
    end
%     fr_range = [fshift(zero_mark:one_thousand_mark),powershift(zero_mark:one_thousand_mark)]'; % 低频信号区
%     fc_range = [fshift(one_thousand_mark:end),powershift(one_thousand_mark:end)]'; % 高频信号区
    [~,fr_index] = max(powershift(zero_mark:one_thousand_mark));
    [~,fc_index] = max(powershift(one_thousand_mark:end));
    fr = fshift(fr_index+zero_mark-1);
    fc = fshift(fc_index+one_thousand_mark-1);
end


function d = find_number(f1,f2)
    % input: 输入低频f1高频f2两个频率值
    % output: 输出对应的键盘号码（0~9）
    fr = [697,770,852,941]';
    fc = [1209 1336 1477]';
    p = NaN;
    q = NaN;
    number = [1, 2, 3;
             4, 5, 6;
             7, 8, 9;
             NaN,0,NaN];
    for i = 1:length(fr)
        if abs(fr(i)-f1)<20
            p = i;
        end
    end
    for j = 1:length(fc)
        if abs(fc(j)-f2)<20
            q = j;
        end
    end
    if (isnan(p) || isnan(q))
        d = NaN;
    else
        d = number(p,q);
    end
end
