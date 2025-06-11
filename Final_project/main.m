clc , clear;
files = dir("Signal*.wav");
for i = 1:length(files)
    [xx, fs] = audioread(files(i).name);
    [~, ~, ~, ~, ~] = toneDetect(xx, fs, 1, files(i).name);
end
